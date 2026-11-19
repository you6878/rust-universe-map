"""Pictures of the Milky Way as the data of the Java program has it, made without running it.

Two sources, both read straight from the data folder:
  - the voxel grids of the galaxy (galaxy/volume/gal-*): light and dust in 3D, with distances;
  - the panorama (tex/skybox/milkyway-*), which the Java program shows while the camera is
    within a kiloparsec of the Sun.

The grids are ray-marched here the way the Java shader does it (emission dimmed by the dust in
front of it); the brightness scale is chosen for the picture, not taken from the shader.

Usage: java_reference.py DATA_DIR OUT_DIR [view ...]
Positions are in parsecs, galactic axes with the Sun at the origin: X to longitude 90, Y to the
galactic north, Z to the galactic centre.
"""
import glob
import gzip
import re
import sys

import numpy as np
from PIL import Image

DATA, OUT = sys.argv[1], sys.argv[2]
SIZE_PC = 26000.0                 # sizePc of the Milky Way object: the unit of the grid bounds
CENTRE = np.array([0.0, 0.0, 8000.0])
DUST_SIGMA = 1000.0               # as in volume.galaxy.fragment.glsl, per unit of grid bounds
W, H, STEPS = 800, 450, 500

# name, where the camera is, what it looks at, horizontal field of view in degrees
VIEWS = [
    ("earth_centre",   (0, 0, 0),              (0, 0, 8000),  90),
    ("earth_l90",      (0, 0, 0),              (8000, 0, 0),  90),
    ("earth_l270",     (0, 0, 0),              (-8000, 0, 0), 90),
    ("above_3kpc",     (0, 3000, 0),           (0, 0, 8000),  90),
    ("half_way",       (0, 300, 4000),         (0, 0, 8000),  90),
    ("edge_outside",   (0, 3000, -22000),      (0, 0, 8000),  60),
    ("top",            (0, 45000, 7999),       (0, 0, 8000),  60),
    ("oblique",        (15000, 15000, -10000), (0, 0, 8000),  60),
]


def pick(pattern):
    """The best quality installed of a path with a '*'."""
    found = glob.glob(f"{DATA}/{pattern}")
    for q in ("high", "med", "low"):
        for f in found:
            if f.endswith(q) or f"-{q}/" in f:
                return f
    return found[0]


def load_grid(folder, kind):
    meta = open(f"{folder}/galaxy_{kind}_metadata.txt").read()
    n = tuple(int(v) for v in re.search(r"^Grid size: (\d+)x(\d+)x(\d+)", meta, re.M).groups())
    b = [float(v) for v in re.findall(r"-?\d+\.\d+", re.search(r"^Volume bounds: (.*)$", meta, re.M).group(1))]
    with gzip.open(f"{folder}/galaxy_{kind}_rgba.raw.gz") as f:
        data = np.frombuffer(f.read(), np.uint8).reshape(n + (4,))
    # Storage axes (W, H, D) are galactic (X to the centre, Z north, -Y): in the axes used here
    # that is (Z, Y, -X).
    return data, np.array(b[:3]) * SIZE_PC, np.array(b[3:]) * SIZE_PC


def rays(pos, target, fov):
    """Unit directions of the pixels of a camera, galactic north up."""
    pos, target = np.array(pos, float), np.array(target, float)
    fwd = (target - pos) / np.linalg.norm(target - pos)
    up = np.array([0.0, 1.0, 0.0]) if abs(fwd[1]) < 0.99 else np.array([0.0, 0.0, 1.0])
    # Longitude grows to the left, as on the sky seen from inside.
    right = np.cross(fwd, up)
    right /= np.linalg.norm(right)
    up = np.cross(right, fwd)
    t = np.tan(np.radians(fov) / 2)
    x = (np.arange(W) + 0.5) / W * 2 - 1
    y = 1 - (np.arange(H) + 0.5) / H * 2
    d = fwd + x[None, :, None] * t * right + y[:, None, None] * t * H / W * up
    return pos, (d / np.linalg.norm(d, axis=-1, keepdims=True)).reshape(-1, 3)


def march(pos, dirs, light, dust, lo, hi):
    """Light of the grid along each ray, dimmed by the dust in front of it."""
    # To the frame of the box: (Z, Y, -X), relative to the galactic centre.
    o = (pos - CENTRE)[[2, 1, 0]] * np.array([1, 1, -1])
    d = dirs[:, [2, 1, 0]] * np.array([1, 1, -1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t0, t1 = (lo - o) / d, (hi - o) / d
    enter = np.maximum(np.minimum(t0, t1).max(axis=1), 0.0)
    leave = np.maximum(t0, t1).min(axis=1)
    length = np.maximum(leave - enter, 0.0)
    step = length / STEPS
    n = np.array(light.shape[:3])
    out = np.zeros((len(d), 3), np.float32)
    through = np.ones(len(d), np.float32)
    jitter = np.random.default_rng(1).random(len(d))
    for k in range(STEPS):
        p = o + d * (enter + (k + jitter) * step)[:, None]
        i = np.clip(((p - lo) / (hi - lo) * n).astype(np.int32), 0, n - 1)
        e = light[i[:, 0], i[:, 1], i[:, 2]].astype(np.float32) / 255.0
        a = dust[i[:, 0], i[:, 1], i[:, 2], 3].astype(np.float32) / 255.0
        dim = np.exp(-DUST_SIGMA * a * (step / SIZE_PC))
        out += e[:, :3] * (e[:, 3] * step / SIZE_PC * through * np.sqrt(dim))[:, None]
        through *= dim
    return out.reshape(H, W, 3), through.reshape(H, W)


def panorama(faces, dirs):
    """The sky box in the directions of the rays. Faces: +X rt, -X lf, +Y up, -Y dn, +Z ft, -Z bk."""
    # The sky box has longitude 90 on its -X side.
    x, y, z = -dirs[:, 0], dirs[:, 1], dirs[:, 2]
    ax, ay, az = np.abs(x), np.abs(y), np.abs(z)
    out = np.zeros((len(dirs), 3), np.uint8)
    # OpenGL cube map rules: for each face, which components run along its u and v.
    rules = [("rt", (ax >= ay) & (ax >= az) & (x > 0), -z, -y, ax), ("lf", (ax >= ay) & (ax >= az) & (x <= 0), z, -y, ax),
             ("up", (ay > ax) & (ay >= az) & (y > 0), x, z, ay), ("dn", (ay > ax) & (ay >= az) & (y <= 0), x, -z, ay),
             ("ft", (az > ax) & (az > ay) & (z > 0), x, -y, az), ("bk", (az > ax) & (az > ay) & (z <= 0), -x, -y, az)]
    for name, mask, u, v, m in rules:
        face = faces[name]
        s = face.shape[0]
        col = np.clip(((u[mask] / m[mask] + 1) / 2 * s).astype(int), 0, s - 1)
        row = np.clip(((v[mask] / m[mask] + 1) / 2 * s).astype(int), 0, s - 1)
        out[mask] = face[row, col]
    return out.reshape(H, W, 3)


def save(name, linear, gain):
    c = linear * gain
    c = (c / (1 + c)) ** (1 / 2.2)
    Image.fromarray((np.clip(c, 0, 1) * 255).astype(np.uint8)).save(f"{OUT}/{name}.png")


folder = pick("default-data/galaxy/volume/gal*")
light, lo, hi = load_grid(folder, "emissive")
dust, _, _ = load_grid(folder, "occlusive")
print(f"grids: {folder}\n  cells {light.shape[:3]}, box {np.round(lo)} to {np.round(hi)} pc from the galactic centre")
sky = pick("default-data/tex/skybox/milkyway*")
faces = {f[-6:-4]: np.asarray(Image.open(f).convert("RGB")) for f in glob.glob(f"{sky}/*.jpg")}
print(f"panorama: {sky}, faces of {faces['ft'].shape[0]} px")

only = sys.argv[3:]
gain = None
for name, pos, target, fov in VIEWS:
    if only and name not in only:
        continue
    pos, dirs = rays(pos, target, fov)
    img, through = march(pos, dirs, light, dust, lo, hi)
    if gain is None:
        gain = 0.6 / max(np.quantile(img, 0.995), 1e-9)       # one exposure for all the views
    save(f"grid_{name}", img, gain)
    line = f"{name}: grid light mean {img.mean():.3e}, darkest dust lets {through.min():.3f} through"
    if np.linalg.norm(pos) < 1000:
        Image.fromarray(panorama(faces, dirs)).save(f"{OUT}/panorama_{name}.png")
        line += "; panorama written (the Java program shows it here)"
    print(line)
