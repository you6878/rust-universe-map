#!/usr/bin/env python3
"""Fits the climate of the planet object to a picture of a real planet.

    fit_climate.py HEIGHTS(.i16.gz without the ending) DAY_CUBEMAP_DIR [--rounds N] [--preview out.jpg]

The planet object (viewer/src/planet.wgsl) computes what grows where from the heights alone:
temperature from latitude and height, rain from the seas upwind and the belts of rising and
sinking air, and from those forest, grass, sand, rock and snow. This script computes the same on
a coarse grid, compares the colors with a cloudless picture of the planet (a cube map of the
Gaia Sky data) over all the land, and searches for the numbers of the model that make them most
alike. It prints the numbers, to be put into the style of the map and the constants of the
shader, and how alike the two are.

The picture is only compared with; nothing of it is kept.
"""
import argparse, glob, gzip, os, sys
import numpy as np
from PIL import Image

W, H = 720, 360


def load_heights(prefix):
    meta = open(prefix + "_metadata.txt").read()
    w, h = [int(v) for v in meta.split("Grid size:")[1].split()[:2]]
    g = np.frombuffer(gzip.open(prefix + ".i16.gz").read(), dtype="<i2").reshape(h, w).astype(np.float32)
    # the mean of the cells of the source in each cell of the coarse grid
    ys = (np.arange(H + 1) * h / H).astype(int); xs = (np.arange(W + 1) * w / W).astype(int)
    out = np.add.reduceat(np.add.reduceat(g, ys[:-1], axis=0), xs[:-1], axis=1)
    return out / (np.diff(ys)[:, None] * np.diff(xs)[None, :])


def load_cube(folder):
    """The cube map as a picture on the grid: pole +Y, longitude 0 along +Z, east towards +X."""
    faces = {}
    for side in ["rt", "lf", "up", "dn", "ft", "bk"]:
        path = glob.glob(os.path.join(folder, f"*_{side}.jpg"))[0]
        faces[side] = np.asarray(Image.open(path).convert("RGB").resize((512, 512)), dtype=np.float32) / 255.0
    lon = (np.arange(W) + 0.5) / W * 2 * np.pi - np.pi
    lat = np.pi / 2 - (np.arange(H) + 0.5) / H * np.pi
    lon, lat = np.meshgrid(lon, lat)
    x, y, z = np.cos(lat) * np.sin(lon), np.sin(lat), np.cos(lat) * np.cos(lon)
    ax, ay, az = np.abs(x), np.abs(y), np.abs(z)
    out = np.zeros((H, W, 3), dtype=np.float32)
    cases = [("rt", (ax >= ay) & (ax >= az) & (x > 0), -z, -y, ax), ("lf", (ax >= ay) & (ax >= az) & (x <= 0), z, -y, ax),
             ("up", (ay > ax) & (ay >= az) & (y > 0), x, z, ay), ("dn", (ay > ax) & (ay >= az) & (y <= 0), x, -z, ay),
             ("ft", (az > ax) & (az > ay) & (z > 0), x, -y, az), ("bk", (az > ax) & (az > ay) & (z <= 0), -x, -y, az)]
    for side, mask, sc, tc, ma in cases:
        u = np.clip(((sc / ma + 1) / 2 * 512).astype(int), 0, 511); v = np.clip(((tc / ma + 1) / 2 * 512).astype(int), 0, 511)
        out[mask] = faces[side][v[mask], u[mask]]
    return out


def smoothstep(a, b, v):
    t = np.clip((v - a) / (b - a), 0, 1)
    return t * t * (3 - 2 * t)


LAT = (np.pi / 2 - (np.arange(H) + 0.5) / H * np.pi)[:, None] * np.ones((1, W))
ROWS, COLS = np.meshgrid(np.arange(H), np.arange(W), indexing="ij")


def shifted(grid, east, north, d):
    """The grid read `d` radians away, `east` and `north` saying which way (their signs may
    differ from row to row)."""
    lat = np.clip(LAT + north * d, -1.55, 1.55)
    rows = np.clip(((np.pi / 2 - lat) / np.pi * H).astype(int), 0, H - 1)
    cols = (COLS + np.rint(east * d / np.cos(LAT).clip(0.05) / (2 * np.pi) * W).astype(int)) % W
    return grid[rows, cols]


# The numbers searched for: name, start, lowest, highest.
PARAMS = [
    ("temp_equator", 32.21, 22, 34), ("temp_pole", -10.52, -40, -10), ("temp_power", 3.089, 1.5, 3.5),
    ("reach", 0.5988, 0.12, 0.6),            # how far upwind the sea still brings rain, radians
    ("near_share", 0.9911, 0.2, 1.0),        # how much less each farther place counts (ratio)
    ("turn_lo", 0.534, 0.3, 0.6), ("turn_hi", 0.8936, 0.5, 0.9),   # latitudes where the wind turns
    ("tilt", 0.05937, 0.0, 1.0),             # how much the trade winds come from the poles' side
    ("monsoon", 0.5035, 0.0, 1.0),           # rain from the sea on the equator's side
    ("rain_lo", 0.1366, 0.0, 0.5), ("rain_hi", 0.9559, 0.4, 1.0),
    ("base", 0.1216, 0.0, 0.3), ("belt", 0.719, 0.0, 0.8), ("sea", 0.8553, 0.2, 1.2), ("sea_in_dry", 0.6659, 0.2, 1.0),
    ("equator_width", 108.1, 10, 120),
    ("cold_bonus", 0.09035, 0.0, 0.6),
    ("plants_lo", 0.1484, 0.05, 0.4), ("plants_hi", 0.5464, 0.25, 0.7), ("trees_lo", 0.4134, 0.3, 0.7), ("trees_hi", 0.9592, 0.45, 1.0),
    ("snow", -13.6, -30, -5),
    ("plants_from", -10.0, -16, 6), ("trees_from", -6.0, -12, 10),   # the coldest that plants and trees stand, C
    ("lapse", 6.5, 4.5, 9.0),
    ("block_lo", 1127, 500, 2500), ("block_hi", 4264, 2500, 5500),
    ("forest_r", 0.1877, 0.02, 0.3), ("forest_g", 0.2391, 0.05, 0.4), ("forest_b", 0.108, 0.01, 0.25),
    ("grass_r", 0.4524, 0.1, 0.6), ("grass_g", 0.2852, 0.1, 0.6), ("grass_b", 0.2182, 0.03, 0.4),
    ("sand_r", 0.717, 0.4, 1.0), ("sand_g", 0.6651, 0.3, 0.9), ("sand_b", 0.4437, 0.15, 0.7),
    ("rock_r", 0.2333, 0.15, 0.7), ("rock_g", 0.2363, 0.15, 0.7), ("rock_b", 0.129, 0.1, 0.6),
]
NAMES = [p[0] for p in PARAMS]


def model(theta, metres):
    """The color of the land everywhere, and how wet and warm it is."""
    p = dict(zip(NAMES, theta))
    alat = np.abs(LAT)
    sign = np.sign(LAT)
    t_sea = p["temp_equator"] - (p["temp_equator"] - p["temp_pole"]) * np.sin(alat) ** p["temp_power"]
    temp = t_sea - p["lapse"] * np.maximum(metres, 0) / 1000.0
    from_east = 1 - smoothstep(p["turn_lo"], p["turn_hi"], alat)
    norm = np.sqrt(1 + p["tilt"] ** 2)
    is_sea = (metres <= 0).astype(np.float32)
    here = np.maximum(metres, 0) * 0.5
    from_sea = np.zeros_like(metres); south_sea = np.zeros_like(metres)
    open_e = np.ones_like(metres); open_w = np.ones_like(metres); open_s = np.ones_like(metres)
    shares = p["near_share"] ** np.arange(4); shares /= shares.sum()
    for k, f in enumerate([0.08, 0.25, 0.55, 1.0]):
        d = p["reach"] * f
        m_e = shifted(metres, 1 / norm, sign * p["tilt"] / norm, d)
        m_w = shifted(metres, -1.0, 0.0, d)
        m_s = shifted(metres, 0.0, -sign, d)
        from_sea += shares[k] * (from_east * open_e * (m_e <= 0) + (1 - from_east) * open_w * (m_w <= 0))
        south_sea += shares[k] * open_s * (m_s <= 0)
        open_e *= 1 - 0.8 * smoothstep(p["block_lo"], p["block_hi"], m_e - here)
        open_w *= 1 - 0.8 * smoothstep(p["block_lo"], p["block_hi"], m_w - here)
        open_s *= 1 - 0.8 * smoothstep(p["block_lo"], p["block_hi"], m_s - here)
    from_sea = np.clip(from_sea + p["monsoon"] * south_sea * from_east, 0, 1)
    rising = 0.5 + 0.5 * np.cos(6 * LAT)
    belts = (rising + (1 - rising) * np.exp(-LAT * LAT * p["equator_width"])) * (1 - 0.6 * smoothstep(1.15, 1.45, alat))
    rain = smoothstep(p["rain_lo"], p["rain_hi"], from_sea)
    wet = np.clip(p["base"] + p["belt"] * belts + p["sea"] * rain * (p["sea_in_dry"] + (1 - p["sea_in_dry"]) * belts), 0, 1)
    enough = wet + p["cold_bonus"] * smoothstep(16, 0, temp)
    plants = smoothstep(p["plants_lo"], p["plants_hi"], enough) * smoothstep(p["plants_from"] - 3, p["plants_from"] + 3, temp)
    trees = smoothstep(p["trees_lo"], p["trees_hi"], enough) * smoothstep(p["trees_from"] - 3, p["trees_from"] + 3, temp)
    col = lambda n: np.array([p[n + "_r"], p[n + "_g"], p[n + "_b"]], dtype=np.float32)
    sand, rock, grass, forest = col("sand"), col("rock"), col("grass"), col("forest")
    w = smoothstep(2, 18, temp)[..., None]
    bare = rock * (1 - w) + sand * w
    g = smoothstep(0.25, 0.5, enough)[..., None]
    grassy = (sand * 0.55 + grass * 0.45) * (1 - g) + grass * g
    wooded = forest * (0.8 + 0.2 * smoothstep(0, 14, temp))[..., None] * (1 - 0.18 * smoothstep(0.7, 0.95, enough))[..., None]
    veg = grassy * (1 - trees[..., None]) + wooded * trees[..., None]
    color = bare * (1 - plants[..., None]) + veg * plants[..., None]
    snow = smoothstep(p["snow"] + 1.5, p["snow"] - 1.5, temp)[..., None]
    color = color * (1 - snow) + np.array([0.95, 0.95, 0.97], dtype=np.float32) * snow
    return color, wet, temp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("heights"); ap.add_argument("cubemap")
    ap.add_argument("--rounds", type=int, default=250); ap.add_argument("--preview")
    args = ap.parse_args()
    metres = load_heights(args.heights)
    photo = load_cube(args.cubemap)
    land = metres > 0
    # The picture and the heights must be the same way round: the seas of the picture are dark blue.
    blue = (photo[..., 2] > photo[..., 0] * 1.3) & (photo.sum(axis=2) < 0.9)
    print(f"the seas of the picture and of the heights agree on {100 * np.mean(blue == ~land):.1f}% of the planet")
    weight = (np.cos(LAT) * land)[..., None]
    weight = weight / weight.sum()

    def loss(theta):
        color, _, _ = model(theta, metres)
        return float((np.abs(color - photo) * weight).sum() / 3)

    from scipy.optimize import differential_evolution
    start = np.array([p[1] for p in PARAMS]); bounds = [(p[2], p[3]) for p in PARAMS]
    print(f"before: the colors of the land differ by {100 * loss(start):.2f}% on average")
    result = differential_evolution(loss, bounds, x0=start, maxiter=args.rounds, popsize=4, tol=1e-7, seed=1, polish=False, updating="deferred", workers=1)
    best = result.x
    color, wet, temp = model(best, metres)
    a, b = (color * np.sqrt(weight * 3)).ravel(), (photo * np.sqrt(weight * 3)).ravel()
    pattern = np.corrcoef(color[land].mean(axis=1), photo[land].mean(axis=1))[0, 1]
    print(f"after:  the colors of the land differ by {100 * loss(best):.2f}% on average; light and dark places agree with r = {pattern:.3f}")
    for name, value in zip(NAMES, best):
        print(f"  {name:14s} {value:.4g}")
    if args.preview:
        sea = np.array([0.03, 0.06, 0.14], dtype=np.float32)
        shown = np.where(land[..., None], color, sea)
        both = np.concatenate([photo, shown], axis=0)
        Image.fromarray((both.clip(0, 1) * 255).astype(np.uint8)).save(args.preview, quality=90)


if __name__ == "__main__":
    main()
