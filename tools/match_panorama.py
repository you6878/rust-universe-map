"""Fills the surroundings of the Sun so that the sky seen from the Earth is the panorama.

The panorama is the picture of the whole sky that the Java program shows from inside the solar
neighbourhood (tex/skybox/milkyway-*). It says how bright and of which color the Milky Way is in
every direction, but not how far anything is. The distances come from elsewhere:
  - the grids of the galaxy as a whole (make_galaxy.py) give the light and the dust far away;
  - the measured dust map (dust_edenhofer.py) gives the dust within 1.25 kpc;
and what is still missing in a direction is put in the box around the Sun: more dust where the
sky is darker than those give, along the measured dust of that direction if there is some and
otherwise in a shell near the edge of the box; light in that shell where the sky is brighter.

So from the Earth the sky is the panorama; from elsewhere, what is added here is where this tool
put it, not where it really is. That is why the additions are a volume of their own, which the
map lets fade away as the camera leaves; the measured dust stays.

Usage: match_panorama.py MAP.json DATA_DIR [PICTURE_OF_THE_FIT.png]
Reads the Milky Way volumes of the map and grids/local_dust_measured.*, writes the additions as
grids/sky_dust.raw.gz, grids/sky_light.raw.gz and grids/sky_metadata.txt next to the map, and
prints the numbers for the map.
"""
import gzip
import json
import os
import re
import sys

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from volumes import LUMA, Field, Grid, light_to_screen, march, panorama, pick_quality, screen_to_light, sky_directions

MAP, DATA = sys.argv[1], sys.argv[2]
BASE = os.path.dirname(os.path.abspath(MAP))
W, H, B_MAX = 1440, 400, 50.0           # the map of the sky the fit is made on: a quarter of a degree
FADE = (34.0, 24.0)                     # nothing is added beyond these latitudes
ROUNDS = 8
MAX_DEPTH = 7.0                         # most optical depth added in a direction
GAMMA = 2.0
INNER = 1800.0                          # parsecs: beyond this only the galaxy as a whole


def smoothstep(a, b, v):
    t = np.clip((v - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


# ---- the map: the galaxy as a whole
the_map = json.load(open(MAP))
volumes = {v["name"]: v for v in the_map["volumes"]}
# The white balance of the camera at the Sun: the viewer multiplies the light by it, so the
# light to make is the panorama's divided by it.
BALANCE = np.array(the_map.get("exposure", {}).get("whiteBalance", [1.0, 1.0, 1.0]), np.float32)
far_light, far_dust = volumes["Milky Way, light"], volumes["Milky Way, dust"]


def grid_of(v, kind):
    return Grid(f"{BASE}/{v[kind]}", f"{BASE}/{v['metadata']}", v["sizePc"], v["centrePc"], v.get("gamma", 1.0))


galaxy_light = grid_of(far_light, "light").light(far_light["lightPerKpc"])
from volumes import load_noise
noise = load_noise(f"{BASE}/{the_map['detail']['noise']}") if "detail" in the_map else None
galaxy_dust = grid_of(far_dust, "dust").dust(far_dust["dustPerKpc"], None if noise is None else (noise, far_dust.get("detailPc", 400.0), far_dust.get("detailStrength", 0.0)))

# ---- the measured dust around the Sun, in optical depth per parsec
meta = open(f"{BASE}/grids/local_dust_measured_metadata.txt").read()
measured_grid = Grid(f"{BASE}/grids/local_dust_measured.raw.gz", f"{BASE}/grids/local_dust_measured_metadata.txt", 1.0, (0, 0, 0), GAMMA)
measured_per_kpc = float(re.search(r"dustPerKpc ([0-9]+\.?[0-9]*)", meta).group(1))
N = measured_grid.n
lo, hi = measured_grid.lo, measured_grid.hi
measured = (measured_grid.data[..., 3].astype(np.float32) / 255.0) ** GAMMA * (measured_per_kpc / 1000.0)

# ---- the sky: what it should be, and the directions
dirs = sky_directions(W, H, B_MAX)
target = screen_to_light(panorama(pick_quality(f"{DATA}/default-data/tex/skybox/milkyway*"), dirs)) / BALANCE
lat = np.degrees(np.arcsin(dirs[:, 1]))
weight = smoothstep(FADE[0], FADE[1], np.abs(lat))
origin = np.zeros(3)

# Light from beyond the box, as it arrives at INNER; and the measured dust in each direction.
beyond, _ = march(origin, dirs, [galaxy_light], [galaxy_dust], start=INNER)
column = np.zeros(len(dirs), np.float32)
t = 0.0
while t < INNER:
    column += Field(measured, lo, hi)(origin + dirs * (t + 3.0)) * 6.0
    t += 6.0
print(f"measured dust: optical depth {np.median(column[np.abs(lat) < 5]):.2f} in the plane (median), {np.median(column[np.abs(lat) > 30]):.2f} above 30 degrees")

# ---- where along a direction the additions go: a shell near the edge of the box
axes = [(np.arange(n, dtype=np.float32) + 0.5) / n * (h - l) + l for n, l, h in zip(N, lo, hi)]


def bump(u, a, b):
    return smoothstep(a, a + 0.15, u) * smoothstep(b, b - 0.13, u)


u = np.linspace(0, 1, 2001)
DUST_SHELL, LIGHT_SHELL = (0.45, 0.98), (0.55, 0.99)
dust_norm, light_norm = np.trapz(bump(u, *DUST_SHELL), u), np.trapz(bump(u, *LIGHT_SHELL), u)


def on_sky(sky_map, lon_deg, lat_deg):
    """Values of a map of the sky (H, W) at directions."""
    col = ((180.0 - lon_deg) % 360.0) / 360.0 * W - 0.5
    row = (B_MAX - lat_deg) / (2 * B_MAX) * H - 0.5
    return map_coordinates(sky_map, [row, col], order=1, mode="nearest")


def fill(extra_depth, extra_light):
    """The box around the Sun: dust (optical depth per parsec) and light (rgb per parsec)."""
    wrap = lambda m: np.concatenate([m[:, -1:], m, m[:, :1]], axis=1)          # longitude goes around
    depth_map, column_map = wrap(extra_depth.reshape(H, W)), wrap(column.reshape(H, W))
    light_maps = [wrap(extra_light[:, k].reshape(H, W)) for k in range(3)]
    dust = np.empty(tuple(N), np.float32)
    light = np.zeros(tuple(N) + (3,), np.float32)
    yy, xx = np.meshgrid(axes[1], axes[2], indexing="ij")
    for i, z in enumerate(axes[0]):                                            # a slab at a time
        d = np.sqrt(z * z + yy * yy + xx * xx)
        lon = np.degrees(np.arctan2(xx, z))
        la = np.degrees(np.arcsin(np.clip(yy / d, -1, 1)))
        col = ((180.0 - lon) % 360.0) / 360.0 * W + 0.5                        # +1 for the wrapped column, -0.5 for the cell centre
        row = np.clip((B_MAX - la) / (2 * B_MAX) * H - 0.5, 0, H - 1)
        look = lambda m: map_coordinates(m, [row, col], order=1, mode="nearest")
        inside = np.abs(la) < B_MAX
        reach = np.minimum(1180.0, 370.0 / np.maximum(np.abs(yy) / d, 1e-3))   # how far the box goes in this direction
        more, have = look(depth_map) * inside, look(column_map)
        along = have / (have + 0.5)                                            # how much of it follows the measured dust
        dust[i] = measured[i] * (1.0 + more * along / np.maximum(have, 1e-3)) + more * (1.0 - along) * bump(d / reach, *DUST_SHELL) / (reach * dust_norm)
        shell = (bump(d / reach, *LIGHT_SHELL) / (reach * light_norm) * inside)[..., None]
        light[i] = np.stack([look(m) for m in light_maps], axis=-1) * shell
    return dust, light


def sky_now(dust, light):
    inner, through = march(origin, dirs, [galaxy_light, Field(light, lo, hi)], [galaxy_dust, Field(dust, lo, hi)], far=INNER, growth=0.015, first=3.0)
    return inner + through * beyond, through


extra_depth = np.zeros(len(dirs), np.float32)
extra_light = np.zeros((len(dirs), 3), np.float32)
dust, light = fill(extra_depth, extra_light)
for round_ in range(ROUNDS + 1):
    now, through = sky_now(dust, light)
    bright_now, bright_target = now @ LUMA, target @ LUMA
    error = np.abs(np.log((bright_now + 0.005) / (bright_target + 0.005)))[weight > 0.99]
    lit = (weight > 0.99) & (bright_target > 0.02)
    tint = lambda v: (v[lit, 0] - v[lit, 2]) / v[lit].sum(axis=1)
    print(f"round {round_}: brightness off by a factor {np.exp(np.median(error)):.2f} (median), {np.exp(np.quantile(error, 0.9)):.2f} (nine in ten); color alike {np.corrcoef(tint(now * BALANCE), tint(target * BALANCE))[0, 1] * 100:.0f} %; "
          f"added depth {extra_depth.mean():.2f}, added light {extra_light.mean():.4f}")
    if round_ == ROUNDS:
        break
    # Too bright: more dust. Too dark: less of the added dust, and where none is left, light.
    extra_depth = np.clip(extra_depth + np.log((bright_now + 1e-4) / (bright_target + 1e-4)), 0.0, MAX_DEPTH) * weight
    # Light is added behind some of the dust of its direction, which dims it: more is asked for there.
    extra_light = np.clip(extra_light + 0.8 * (target - now) / np.clip(np.sqrt(through), 0.3, 1.0), 0.0, None) * weight[:, None]
    dust, light = fill(extra_depth, extra_light)


def write(name, density, rgb):
    top = float(np.quantile(density[density > 0], 0.9998)) if (density > 0).any() else 1.0
    out = np.empty(tuple(N) + (4,), np.uint8)
    out[..., :3] = rgb
    out[..., 3] = (np.clip(density / top, 0.0, 1.0) ** (1.0 / GAMMA) * 255.0 + 0.5).astype(np.uint8)
    with gzip.open(f"{BASE}/grids/{name}.raw.gz", "wb", compresslevel=4) as f:
        f.write(out.tobytes())
    return top * 1000.0


dust_per_kpc = write("sky_dust", np.maximum(dust - measured, 0.0), np.array([140, 97, 61], np.uint8))
strength = light.max(axis=-1)
light_per_kpc = write("sky_light", strength, (light / np.maximum(strength, 1e-12)[..., None] * 255.0 + 0.5).astype(np.uint8))
with open(f"{BASE}/grids/sky_metadata.txt", "w") as f:
    f.write("# Dust and light around the Sun that make the sky of the Earth the Milky Way panorama of Gaia Sky,\n")
    f.write("# on top of the galaxy as a whole and of the measured dust. Made by tools/match_panorama.py.\n")
    f.write("# Axes: (galactic centre, north, longitude 90).\n")
    f.write(f"Grid size: {N[0]}x{N[1]}x{N[2]}\n")
    f.write(f"Volume bounds: [{lo[0]}, {lo[1]}, {lo[2]}] to [{hi[0]}, {hi[1]}, {hi[2]}]\n")
    f.write(f"# Bounds are in parsecs: sizePc 1, centrePc [0, 0, 0], gamma {GAMMA:.0f}, lightPerKpc {light_per_kpc:.4g}, dustPerKpc {dust_per_kpc:.4g}.\n")
print(f'for the map: "sizePc": 1, "centrePc": [0, 0, 0], "gamma": {GAMMA:.0f}, "lightPerKpc": {light_per_kpc:.4g}, "dustPerKpc": {dust_per_kpc:.4g}')
if len(sys.argv) > 3:                   # a picture of the fit: the panorama above, the result below
    Image.fromarray(np.vstack([light_to_screen(target * BALANCE).reshape(H, W, 3), light_to_screen(now * BALANCE).reshape(H, W, 3)])).save(sys.argv[3])
