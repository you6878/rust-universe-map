"""What the tools share: reading the voxel grids of the map, marching rays through them the way
the viewer's shader does, and reading the panorama of the Java program.

Axes everywhere: galactic, parsecs, the Sun at the origin; X to longitude 90, Y to the galactic
north, Z to the galactic centre.
"""
import glob
import gzip
import re

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates

REDDENING = np.array([0.95, 1.0, 1.05], np.float32)     # as in volume.wgsl


class Grid:
    """A voxel grid of the map: RGBA bytes on a box."""

    def __init__(self, raw, metadata, size_pc, centre_pc, gamma=2.0):
        meta = open(metadata).read()
        self.n = tuple(int(v) for v in re.search(r"^Grid size: (\d+)x(\d+)x(\d+)", meta, re.M).groups())
        b = [float(v) for v in re.findall(r"-?\d+\.?\d*", re.search(r"^Volume bounds: (.*)$", meta, re.M).group(1))]
        # The box has its axes in the order (Z, Y, X) of the galactic ones.
        centre = np.array(centre_pc, float)[::-1]
        self.lo, self.hi = np.array(b[:3]) * size_pc + centre, np.array(b[3:]) * size_pc + centre
        with gzip.open(raw) as f:
            self.data = np.frombuffer(f.read(), np.uint8).reshape(self.n + (4,))
        self.gamma = gamma
        self.cell = float(((self.hi - self.lo) / self.n).min())

    def sample(self, p, channels):
        """The bytes at points p (galactic x, y, z), smoothed between cells as the GPU does; 0..1."""
        c = (p[:, ::-1] - self.lo) / (self.hi - self.lo) * self.n - 0.5
        inside = np.all((c > -0.5) & (c < np.array(self.n) - 0.5), axis=1)
        return [map_coordinates(self.data[..., ch], c.T, order=1, mode="nearest", output=np.float32) / 255.0 * inside for ch in channels]

    def light(self, per_kpc):
        """Light per parsec at points, for march(): rgb."""
        def f(p):
            r, g, b, a = self.sample(p, (0, 1, 2, 3))
            return np.stack([r, g, b], axis=1) * (a ** self.gamma * per_kpc / 1000.0)[:, None]
        return f

    def dust(self, per_kpc, detail=None):
        """Optical depth per parsec at points, for march(). `detail` = (noise cube, parsecs it
        repeats after, strength) adds the fine clumps the way the shader does."""
        if detail is None or detail[2] <= 0:
            return lambda p: self.sample(p, (3,))[0] ** self.gamma * per_kpc / 1000.0
        cube, period, strength = detail

        def f(p):
            q = p[:, ::-1] / period                               # the box order of the axes, as in the shader
            two = noise_at(cube, q) + noise_at(cube, q[:, [1, 2, 0]] * 0.2967 + 0.37)
            bell = (two - 256.0 / 255.0) * (255.0 / 48.0) * 0.70711
            return self.sample(p, (3,))[0] ** self.gamma * per_kpc / 1000.0 * np.exp(strength * bell - 0.5 * strength * strength)
        return f


def load_noise(path):
    """The cube of noise of the map: bytes, as many cells along each side."""
    with gzip.open(path) as f:
        data = np.frombuffer(f.read(), np.uint8)
    side = round(len(data) ** (1 / 3))
    return data.reshape(side, side, side)


def noise_at(cube, q):
    """The noise at texture coordinates q (x, y, z), 0..1, repeating: as a GPU reads a texture
    whose first axis in the file is its depth."""
    side = cube.shape[0]
    index = (q[:, ::-1] * side - 0.5).T
    return map_coordinates(cube, index, order=1, mode="grid-wrap", output=np.float32) / 255.0


class Field:
    """Numbers on the cells of a box, not yet turned into bytes: one value or a color per cell, per parsec."""

    def __init__(self, values, lo, hi):
        self.values, self.lo, self.hi = values, np.array(lo, float), np.array(hi, float)
        self.n = np.array(values.shape[:3])

    def __call__(self, p):
        c = (p[:, ::-1] - self.lo) / (self.hi - self.lo) * self.n - 0.5
        inside = np.all((c > -0.5) & (c < self.n - 0.5), axis=1)
        if self.values.ndim == 3:
            return map_coordinates(self.values, c.T, order=1, mode="nearest") * inside
        return np.stack([map_coordinates(self.values[..., k], c.T, order=1, mode="nearest") for k in range(3)], axis=1) * inside[:, None]


def march(origin, dirs, lights, dusts, start=0.0, far=45000.0, growth=0.03, first=4.0):
    """Light along rays between two distances, as the viewer adds it up.

    lights and dusts are functions of points that give the light (rgb) and the optical depth per
    parsec. Returns (light rgb, what the dust lets through in rgb), per ray. Steps grow with the
    distance, as in the shader.
    """
    n = len(dirs)
    out = np.zeros((n, 3), np.float32)
    through = np.ones((n, 3), np.float32)
    t = start
    jitter = np.random.default_rng(3).random(n).astype(np.float32)
    while t < far:
        step = min(t * growth + first, far - t)
        p = origin + dirs * (t + step * jitter)[:, None]
        shine = np.zeros((n, 3), np.float32)
        depth = np.zeros(n, np.float32)
        for f in lights:
            shine += f(p)
        for f in dusts:
            depth += f(p)
        dim = np.exp(-depth[:, None] * step * REDDENING)
        out += shine * step * through * np.sqrt(dim)
        through *= dim
        t += step
    return out, through


def sky_directions(width, height, b_max):
    """Directions of a map of the sky: longitude 180 .. -180 left to right, latitude b_max .. -b_max."""
    lon = np.radians(180.0 - (np.arange(width) + 0.5) / width * 360.0)
    lat = np.radians(b_max - (np.arange(height) + 0.5) / height * 2 * b_max)
    lon, lat = np.meshgrid(lon, lat)
    return np.stack([np.cos(lat) * np.sin(lon), np.sin(lat), np.cos(lat) * np.cos(lon)], axis=-1).reshape(-1, 3)


def panorama(folder, dirs):
    """The sky box of the Java program in the given directions, as bytes."""
    faces = {f[-6:-4]: np.asarray(Image.open(f).convert("RGB")) for f in glob.glob(f"{folder}/*.jpg")}
    # The sky box has longitude 90 on its -X side.
    x, y, z = -dirs[:, 0], dirs[:, 1], dirs[:, 2]
    ax, ay, az = np.abs(x), np.abs(y), np.abs(z)
    out = np.zeros((len(dirs), 3), np.uint8)
    # OpenGL cube map rules: +X rt, -X lf, +Y up, -Y dn, +Z ft, -Z bk.
    rules = [("rt", (ax >= ay) & (ax >= az) & (x > 0), -z, -y, ax), ("lf", (ax >= ay) & (ax >= az) & (x <= 0), z, -y, ax),
             ("up", (ay > ax) & (ay >= az) & (y > 0), x, z, ay), ("dn", (ay > ax) & (ay >= az) & (y <= 0), x, -z, ay),
             ("ft", (az > ax) & (az > ay) & (z > 0), x, -y, az), ("bk", (az > ax) & (az > ay) & (z <= 0), -x, -y, az)]
    for name, mask, u, v, m in rules:
        face = faces[name]
        s = face.shape[0]
        out[mask] = face[np.clip(((v[mask] / m[mask] + 1) / 2 * s).astype(int), 0, s - 1), np.clip(((u[mask] / m[mask] + 1) / 2 * s).astype(int), 0, s - 1)]
    return out


LUMA = np.array([0.2126, 0.7152, 0.0722], np.float32)          # as in volume.wgsl


def screen_to_light(pixels):
    """From the bytes of a picture to the light the viewer would need for them: the reverse of
    its tone curve (sRGB of (l / (1 + brightness of l)) ** 1.35)."""
    v = pixels.astype(np.float32) / 255.0
    lin = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
    c = lin ** (1 / 1.35)
    return c / (1.0 - np.clip(c @ LUMA, 0.0, 0.985))[..., None]


def light_to_screen(light):
    c = np.clip(light / (1.0 + light @ LUMA)[..., None], 0.0, 1.0) ** 1.35
    v = np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)
    return (np.clip(v, 0, 1) * 255 + 0.5).astype(np.uint8)


def pick_quality(pattern):
    """The best quality installed of a path with a '*' (-high, -med, -low)."""
    found = glob.glob(pattern)
    for q in ("high", "med", "low"):
        for f in found:
            if f.endswith(q):
                return f
    return found[0]
