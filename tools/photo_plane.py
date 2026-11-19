"""Makes the face-on map of a galaxy's disc from its photograph, finer than its voxel grid.

The light grid of tools/photo_galaxy.py has cells of a few hundred parsecs, and the dust lanes in
it are filled in (they are put in the dust grid). The star formula (viewer/src/procedural.rs,
GridGalaxy) reads this map besides the grid: within each cell of the grid it spreads the stars as
the photograph does, lanes and all, so that seen from where the photograph was taken the stars
look like it.

The disc is found and stretched into a circle exactly as photo_galaxy.py does it (same ellipse,
same axes: first the long axis, second the short one), so the map lies on the grid's plane.

With --map (the map's universe.json) and --volume (the name of the galaxy's light volume in it)
the map is then fitted: the stars of the formula (the grid's cells, each spread by the map) are
seen from the Earth on the photograph's pixels, and the map is corrected where they are too bright
or too faint, a few times over. The photograph is what the galaxy looks like from here, height and
tilt and all, so this takes in what a face-on map alone cannot (the thickness of the disc, the
round bulge that the stretch draws out).

Usage: photo_plane.py PICTURE OUT_DIR NAME --radius-kpc R [--n 1280] [--map FILE --volume NAME]
Writes NAME_plane.png (16-bit, light in linear units, the brightest 1) and NAME_plane_metadata.txt.
"""
import gzip
import json
import argparse

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates

args = argparse.ArgumentParser()
args.add_argument("picture")
args.add_argument("out")
args.add_argument("name")
args.add_argument("--radius-kpc", type=float, required=True, help="as given to photo_galaxy.py")
args.add_argument("--n", type=int, default=1280, help="pixels across the map")
args.add_argument("--map", help="universe.json, to fit the map to the photograph")
args.add_argument("--volume", help="the galaxy's light volume in it")
args.add_argument("--rounds", type=int, default=12)
a = args.parse_args()

HALF = a.radius_kpc * 1.2                           # kpc: the grid's plane reaches as far
N = a.n
CELL = 2 * HALF / N


def axis(n, half):
    return (np.arange(n, dtype=np.float32) + 0.5) / n * 2 * half - half


def smoothstep(lo, hi, v):
    t = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3 - 2 * t)


# ---- the ellipse in the photograph (as photo_galaxy.py)
photo = (np.asarray(Image.open(a.picture).convert("RGB"), np.float32) / 255.0) ** 2.2
lum = photo.mean(axis=-1)
smooth = gaussian_filter(lum, 6)
cy, cx = np.unravel_index(smooth.argmax(), smooth.shape)
ys, xs = np.nonzero(smooth > 0.02 * smooth.max())
cov = np.cov(np.stack([xs - cx, ys - cy]))
values, vectors = np.linalg.eigh(cov)
long_axis = vectors[:, 1] * np.sign(vectors[0, 1] or 1.0)
short_axis = np.array([-long_axis[1], long_axis[0]])
along = (xs - cx) * long_axis[0] + (ys - cy) * long_axis[1]
reach_px, flat = np.quantile(np.abs(along), 0.995), float(np.sqrt(values[0] / values[1]))
kpc_per_px = a.radius_kpc / reach_px

# ---- the disc as from above, on the grid's plane, without the lanes filled in
u = axis(N, HALF)[:, None]
v = axis(N, HALF)[None, :]
col = cx + (u * long_axis[0] + v * flat * short_axis[0]) / kpc_per_px
row = cy + (u * long_axis[1] + v * flat * short_axis[1]) / kpc_per_px
plane = map_coordinates(lum, [row, col], order=1, cval=0.0) * flat
plane *= smoothstep(HALF, HALF * 0.9, np.sqrt(u * u + v * v))
plane = np.clip(plane / plane.max(), 0.0, 1.0)


def rot(which, deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return {"y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[which]


if a.map:
    import os
    u_map = json.load(open(a.map))
    base = os.path.dirname(a.map)
    vol = next(v for v in u_map["volumes"] if v["name"] == a.volume)
    meta = open(os.path.join(base, vol["metadata"])).read()
    n0, n1, n2 = [int(x) for x in meta.split("Grid size:")[1].split("\n")[0].strip().split("x")]
    raw = np.frombuffer(gzip.open(os.path.join(base, vol["light"])).read(), np.uint8).reshape(n0, n1, n2, 4)
    dens = (raw[..., 3] / 255.0) ** vol["gamma"]
    # the grid's stars, as many points: in their cells at random, each with the map's factor
    rng = np.random.default_rng(5)
    P_N = 6_000_000
    cells = rng.choice(dens.size, P_N, p=(dens / dens.sum()).ravel())
    i0, rest = np.divmod(cells, n1 * n2)
    i1, i2 = np.divmod(rest, n2)
    half_h = float(meta.split("Volume bounds:")[1].split("]")[0].split(",")[1])
    q = np.stack([(i0 + rng.random(P_N)) / n0 * 2 * HALF - HALF, (i1 + rng.random(P_N)) / n1 * 2 * half_h - half_h, (i2 + rng.random(P_N)) / n2 * 2 * HALF - HALF], -1)
    pix = np.clip(((q[:, 0] + HALF) / (2 * HALF) * N).astype(int), 0, N - 1) * N + np.clip(((q[:, 2] + HALF) / (2 * HALF) * N).astype(int), 0, N - 1)
    # seen from the Earth (the Sun), on the photograph's pixels
    to_eq = rot("y", 90 + 192.85948) @ rot("z", 90 - 27.12825) @ rot("y", -32.93192)
    centre = to_eq @ np.array(vol["centrePc"])
    sight = centre / np.linalg.norm(centre)
    north = np.array([0.0, 1.0, 0.0]) - sight[1] * sight
    north /= np.linalg.norm(north)
    east = np.cross(north, sight)
    e = (np.array(vol["centrePc"]) + (q * 1000.0) @ np.array(vol["axes"])) @ to_eq.T
    d = e / np.linalg.norm(e, axis=1, keepdims=True)
    rad_per_px = kpc_per_px * 1000.0 / np.linalg.norm(centre)
    col = np.round(cx - (d @ east) / rad_per_px).astype(int)
    row = np.round(cy - (d @ north) / rad_per_px).astype(int)
    H, W = lum.shape
    seen = (col >= 0) & (col < W) & (row >= 0) & (row < H)
    at = np.where(seen, row * W + col, 0)
    target = gaussian_filter(lum, 1.0).ravel()
    # the factor of the map, fitted: each point carries the factor of its pixel of the map
    factor = np.ones(N * N)
    count = np.bincount(pix, minlength=N * N)
    for k in range(a.rounds):
        w = factor[pix] * seen
        F = gaussian_filter(np.bincount(at, weights=w, minlength=H * W).reshape(H, W), 1.0).ravel()
        F *= target.sum() / F.sum()
        want = np.where(F > 1e-9, target / np.maximum(F, 1e-9), 1.0)
        # each pixel of the map: the mean correction of its points, half of it, smoothed a little
        fix = np.bincount(pix, weights=want[at] * seen + (1 - seen), minlength=N * N) / np.maximum(count, 1)
        fix = gaussian_filter(np.where(count > 0, fix, 1.0).reshape(N, N), 0.7).ravel()
        factor = np.clip(factor * fix ** 0.7, 0.0, 30.0)
        ok = target > 0.02 * target.max()
        corr = np.corrcoef(np.sqrt(F[ok]), np.sqrt(target[ok]))[0, 1]
        print(f"round {k + 1}: like the photograph {corr:.3f}")
    # the map: the grid seen face on, times the factor (the formula takes the factor back out)
    column = dens.sum(axis=1)
    grid_face = map_coordinates(column, [((u + HALF) / (2 * HALF) * n0 - 0.5) + 0 * v, ((v + HALF) / (2 * HALF) * n2 - 0.5) + 0 * u], order=1, mode="nearest")
    plane = grid_face * factor.reshape(N, N)
    plane = np.clip(plane / plane.max(), 0.0, 1.0)
Image.fromarray((plane * 65535 + 0.5).astype(np.uint16)).save(f"{a.out}/{a.name}_plane.png")
with open(f"{a.out}/{a.name}_plane_metadata.txt", "w") as f:
    f.write(f"# Made by tools/photo_plane.py from {a.picture.split('/')[-1]}: the disc as from above, rows along its long axis, columns along its short one.\n")
    f.write(f"Plane size: {N}x{N}\n")
    f.write(f"Half width kpc: {HALF}\n")
    f.write(f"# Where the photograph has the galaxy (for comparing with it): its middle and its scale.\n")
    f.write(f"Photograph middle px: {cx} {cy}\n")
    f.write(f"Photograph kpc per px: {kpc_per_px}\n")
print(f"wrote {a.out}/{a.name}_plane.png: {N}x{N}, {CELL * 1000:.0f} pc a pixel; the photograph has {kpc_per_px * 1000:.0f} pc a pixel")
