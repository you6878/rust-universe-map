"""Makes the voxel grids of a spiral galaxy from a photograph of it as it is seen in the sky.

A disc seen at a slant is an ellipse. The tool finds the ellipse in the photograph (its middle,
the direction of its long axis and how flat it looks), stretches it back into a circle to get
the galaxy as from above, and gives that a thickness the way make_galaxy.py does for the Milky
Way: a smooth haze of old stars, what stands out of it in a thin layer, a round bulge, and dust
where the photograph has dark lanes. How flat the ellipse looks also says how much the disc is
tilted, and with the place of the galaxy in the sky that gives where its three axes point.

What the photograph cannot say: which things are in front of which. The far half and the near
half are told apart only by `--near`, and the bulge, which is round and not flat, is made round
again from how bright the middle is.

Usage: photo_galaxy.py PICTURE OUT_DIR NAME --radius-kpc R --ra DEG --dec DEG --distance-pc D
                       [--near nw|ne|sw|se] [--bulge-gain G] [--preview DIR]
The photograph must have north up and east to the left. Writes NAME_light.raw.gz,
NAME_dust.raw.gz and NAME_metadata.txt, and prints the volume to put in the map.
"""
import argparse
import gzip
import json

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates

args = argparse.ArgumentParser()
args.add_argument("picture")
args.add_argument("out")
args.add_argument("name")
args.add_argument("--radius-kpc", type=float, required=True, help="how far the disc of the photograph reaches along its long axis")
args.add_argument("--ra", type=float, required=True)
args.add_argument("--dec", type=float, required=True)
args.add_argument("--distance-pc", type=float, required=True)
args.add_argument("--near", default="nw", help="the side of the short axis that is nearer to us")
args.add_argument("--bulge-gain", type=float, default=1.0)
args.add_argument("--preview")
a = args.parse_args()

N_PLANE, LIGHT_N, DUST_N = 320, 48, 16
HALF = a.radius_kpc * 1.2                           # kpc, in the plane
LIGHT_HALF, DUST_HALF = 3.0, 0.5                    # kpc, above and below the plane
GAMMA = 2.0
CELL = 2 * HALF / N_PLANE


def axis(n, half):
    return (np.arange(n, dtype=np.float32) + 0.5) / n * 2 * half - half


def smoothstep(lo, hi, v):
    t = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def wide(v, kpc):
    return gaussian_filter(v, kpc / CELL)


# ---- the ellipse in the photograph
photo = (np.asarray(Image.open(a.picture).convert("RGB"), np.float32) / 255.0) ** 2.2
lum = photo.mean(axis=-1)
smooth = gaussian_filter(lum, 6)
cy, cx = np.unravel_index(smooth.argmax(), smooth.shape)
# The outline of the disc: everything brighter than a little of the sky's darkness.
ys, xs = np.nonzero(smooth > 0.02 * smooth.max())
cov = np.cov(np.stack([xs - cx, ys - cy]))
values, vectors = np.linalg.eigh(cov)
long_axis = vectors[:, 1] * np.sign(vectors[0, 1] or 1.0)        # in the picture: (to the right, down)
short_axis = np.array([-long_axis[1], long_axis[0]])
along, across = (xs - cx) * long_axis[0] + (ys - cy) * long_axis[1], (xs - cx) * short_axis[0] + (ys - cy) * short_axis[1]
# How flat: an evenly filled ellipse is spread along its two axes in the ratio of their lengths.
reach_px, flat = np.quantile(np.abs(along), 0.995), float(np.sqrt(values[0] / values[1]))
kpc_per_px = a.radius_kpc / reach_px
tilt = np.degrees(np.arccos(flat))
angle = np.degrees(np.arctan2(-long_axis[0], -long_axis[1])) % 180       # of the long axis, from north through east
print(f"the disc in the photograph: middle at ({cx}, {cy}) px, long axis {reach_px:.0f} px = {a.radius_kpc} kpc, "
      f"{angle:.0f} degrees from north through east, {flat:.2f} times as wide as long: tilted by {tilt:.0f} degrees")

# ---- the galaxy as from above: the ellipse stretched into a circle. First axis of the plane:
# the long axis; third: the short one. A tilted disc looks brighter, as more of it is behind
# each point: that is taken out again.
u = axis(N_PLANE, HALF)[:, None]
v = axis(N_PLANE, HALF)[None, :]
col = cx + (u * long_axis[0] + v * flat * short_axis[0]) / kpc_per_px
row = cy + (u * long_axis[1] + v * flat * short_axis[1]) / kpc_per_px
blur = max(CELL / kpc_per_px / 2.5, 0.5)
plane = np.stack([map_coordinates(gaussian_filter(photo[..., c], blur), [row, col], order=1, cval=0.0) for c in range(3)], axis=-1) * flat
r = np.sqrt(u * u + v * v)
plane *= smoothstep(HALF, HALF * 0.9, r)[..., None]
brightness = plane.mean(axis=-1)

# The bulge is round, not flat, so the stretch has drawn it out along the short axis. It is put
# back as it is along the long axis, where the tilt does not show, the same all around.
profile_r = np.arange(0, N_PLANE // 2) * CELL
on_axis = np.abs(v) < 2 * CELL + 0 * u
profile = np.stack([[np.median(plane[..., c][on_axis & (np.abs(np.abs(u) - d) < CELL)]) for d in profile_r] for c in range(3)], axis=-1)
# It ends where the light has come down to the level of the disc: a few times what the disc
# has half-way out, and never more than a quarter of the way out.
along_axis = np.nan_to_num(profile.mean(axis=-1))
disc_level = np.median(along_axis[(profile_r > 0.35 * a.radius_kpc) & (profile_r < 0.7 * a.radius_kpc)])
bulge_r = min(profile_r[np.argmax(along_axis < 2.5 * disc_level)], 0.25 * a.radius_kpc)
print("light along the long axis, every 2 kpc:", " ".join(f"{along_axis[int(d / CELL)]:.3f}" for d in np.arange(0, a.radius_kpc, 2.0)))
round_middle = np.stack([np.interp(r, profile_r, profile[:, c]) for c in range(3)], axis=-1)
inside = smoothstep(2.2 * bulge_r, 1.2 * bulge_r, r)[..., None]
plane = plane * (1 - inside) + round_middle * inside
bulge = wide(smoothstep(1.6 * bulge_r, 0.8 * bulge_r, r), 0.3)
lum = plane.mean(axis=-1)
print(f"the bulge reaches {bulge_r:.1f} kpc from the middle")

# Dust lanes: darker than their surroundings.
around = wide(lum, 0.6)
# (Only clearly darker counts: a photograph has grain.)
lanes = wide(np.clip((1.0 - lum / (around + 1e-5) - 0.1) * 4.0, 0.0, 1.0), 0.15) * (1.0 - bulge) * smoothstep(0.03, 0.1, around / around.max())
keep = 1.0 - smoothstep(0.1, 0.5, lanes)
filled = np.stack([gaussian_filter(plane[..., c] * keep, 3.0) / (gaussian_filter(keep, 3.0) + 1e-3) for c in range(3)], axis=-1)
shine = plane * keep[..., None] + filled * (1.0 - keep[..., None])
haze = np.minimum(shine, np.stack([wide(shine[..., c], 0.9) for c in range(3)], axis=-1))
arms = shine - haze
haze_disc, arms_disc, bulge_light = haze * (1 - bulge[..., None]), arms * (1 - bulge[..., None]), shine * bulge[..., None]


def layer(n, half, profile):
    """A vertical profile on the cells of a grid, with an integral of 1 over the height (in kpc)."""
    y = axis(n, half)
    p = profile(np.abs(y)) * smoothstep(half, half * 0.8, np.abs(y))
    return (p / (p.sum() * 2 * half / n)).astype(np.float32)


def write(kind, density, rgb, top, n, half):
    out = np.empty(density.shape + (4,), np.uint8)
    out[..., :3] = (np.clip(rgb, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    out[..., 3] = (np.clip(density / top, 0.0, 1.0) ** (1.0 / GAMMA) * 255.0 + 0.5).astype(np.uint8)
    with gzip.open(f"{a.out}/{a.name}_{kind}.raw.gz", "wb", compresslevel=4) as f:
        f.write(out.tobytes())
    with open(f"{a.out}/{a.name}_{kind}_metadata.txt", "w") as f:
        f.write(f"# Made by tools/photo_galaxy.py from {a.picture.split('/')[-1]}. Axes: (the long axis of the disc, its pole, its short axis).\n")
        f.write(f"Grid size: {N_PLANE}x{n}x{N_PLANE}\n")
        f.write(f"Volume bounds: [{-HALF}, {-half}, {-HALF}] to [{HALF}, {half}, {HALF}]\n")
        f.write(f"# Bounds are in kpc: sizePc 1000, gamma {GAMMA:.0f}.\n")


# ---- light: each layer of the plane times how it is spread in height
thin = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / 0.12))
thick = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / 0.35))
tall = max(bulge_r * 0.7, 0.4)
round_ = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / (0.45 * tall)))
reach = wide(bulge, 0.5)
reach = np.maximum(min(2.2 * tall, LIGHT_HALF * 0.9) * reach / reach.max(), 1e-3)
rounded = smoothstep(reach[:, None, :], 0.45 * reach[:, None, :], np.abs(axis(LIGHT_N, LIGHT_HALF))[None, :, None] + 0 * reach[:, None, :])
light = (haze_disc[:, None, :, :] * thick[None, :, None, None] + arms_disc[:, None, :, :] * thin[None, :, None, None]
         + a.bulge_gain * bulge_light[:, None, :, :] * round_[None, :, None, None] * rounded[..., None])
strength = light.max(axis=-1)
light_top = float(np.quantile(strength[strength > 0], 0.9995))
write("light", strength, light / np.maximum(strength, 1e-9)[..., None], light_top, LIGHT_N, LIGHT_HALF)
del light, strength

# ---- dust: spread out through the disc, and thick where the photograph has dark lanes
rng = np.random.default_rng(abs(hash(a.name)) % (2 ** 32))
shape = (N_PLANE, DUST_N, N_PLANE)
white = np.fft.rfftn(rng.standard_normal(shape, dtype=np.float32))
k = np.sqrt(sum(np.meshgrid(np.fft.fftfreq(shape[0]) ** 2, (np.fft.fftfreq(shape[1]) * N_PLANE / DUST_N * DUST_HALF / HALF * 0.5) ** 2, np.fft.rfftfreq(shape[2]) ** 2, indexing="ij")))
k[0, 0, 0] = 1.0
noise = np.fft.irfftn(white * k ** (-1.4), s=shape).astype(np.float32)
noise = (noise - noise.mean()) / noise.std()
disc = smoothstep(0.03, 0.12, wide(lum, 1.0) / lum.max()) * (0.4 + 0.6 * smoothstep(0.5 * bulge_r, 1.5 * bulge_r, r))
in_plane = 0.7 * disc + 10.0 * lanes
height = layer(DUST_N, DUST_HALF, lambda y: np.exp(-y / 0.1)) * 0.2
dust = in_plane[:, None, :] * height[None, :, None] * np.exp(0.8 * noise - 0.32)
dust_top = float(np.quantile(dust[dust > 0], 0.9995))
write("dust", dust, np.array([0.55, 0.38, 0.24], np.float32), dust_top, DUST_N, DUST_HALF)

# ---- where the galaxy is and how it lies, in the axes of the map (galactic, parsecs)
def rot(which, deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return {"y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[which]


# From the galactic axes of the map to the equatorial ones (+Y north pole, +Z equinox, +X RA 90).
to_equatorial = rot("y", 90 + 192.85948) @ rot("z", 90 - 27.12825) @ rot("y", -32.93192)
ra, dec = np.radians(a.ra), np.radians(a.dec)
sight = np.array([np.cos(dec) * np.sin(ra), np.sin(dec), np.cos(dec) * np.cos(ra)])
north = np.array([0.0, 1.0, 0.0]) - sight[1] * sight
north /= np.linalg.norm(north)
east = np.cross(north, sight)
# The picture has north up and east to the left.
sky = lambda picture_vector: -picture_vector[0] * east - picture_vector[1] * north
long_sky, short_sky = sky(long_axis), sky(short_axis)
# The short axis leans towards us on one side and away on the other.
near = {"n": north, "s": -north}[a.near[0]] + {"e": east, "w": -east}[a.near[1]]
leans_away = 1.0 if short_sky @ near < 0 else -1.0
short_in_space = flat * short_sky + leans_away * np.sqrt(1 - flat * flat) * sight
pole = np.cross(long_sky, short_in_space)
axes = [(to_equatorial.T @ w).round(6).tolist() for w in (long_sky, pole / np.linalg.norm(pole), short_in_space)]
centre = (to_equatorial.T @ (sight * a.distance_pc)).round(1).tolist()

# How bright: as much light from above as the disc of the Milky Way of the map gives.
disc_light = np.median(haze_disc.mean(axis=-1)[(r > 0.3 * a.radius_kpc) & (r < 0.7 * a.radius_kpc)])
volume = {
    "name": a.name,
    "note": f"Made by tools/photo_galaxy.py from the photograph {a.picture.split('/')[-1]}: its ellipse stretched into a disc, tilted by {tilt:.0f} degrees. Which half is nearer and the roundness of the bulge are not in a photograph.",
    "light": f"grids/{a.name}_light.raw.gz", "metadata": f"grids/{a.name}_light_metadata.txt",
    "sizePc": 1000, "centrePc": centre, "axes": axes, "gamma": 2, "lightPerKpc": round(light_top * 0.009 / max(disc_light, 1e-9), 4),
}
dust_volume = {
    "name": f"{a.name}, dust", "dust": f"grids/{a.name}_dust.raw.gz", "metadata": f"grids/{a.name}_dust_metadata.txt",
    "sizePc": 1000, "centrePc": centre, "axes": axes, "gamma": 2, "dustPerKpc": round(dust_top, 2), "dustGlowPerKpc": 0.02, "detailPc": 400, "detailStrength": 1.0,
}
print("for the map, in \"volumes\":")
print(json.dumps([volume, dust_volume]))

if a.preview:
    for name, picture in (("above", shine / shine.max()), ("lanes", lanes), ("bulge", bulge)):
        Image.fromarray((np.clip(picture, 0, 1) ** (1 / 2.2) * 255).astype(np.uint8)).save(f"{a.preview}/{a.name}_{name}.png")
