"""Makes the voxel grids of the Milky Way as a whole from a picture of it seen from above.

The picture is the face-on map by R. Hurt (NASA/JPL-Caltech/ESO), an artist's impression drawn
from the measured positions of the bar, the arms and the Sun; there is no photograph of the
galaxy from outside. https://www.eso.org/public/images/eso1339g/ (5600 x 5600 px).

The picture gives where things are in the plane of the galaxy: its light (the haze of the disc,
the star clouds of the arms, the pink gas clouds, the yellow bar) and its dust (the brown lanes).
How thick each of them is above and below the plane is a model. The dust leaves out the
surroundings of the Sun, which the measured dust map fills (see dust_edenhofer.py).

Usage: make_galaxy.py PICTURE OUT_DIR [TUNING_DIR]
Writes mw_light.raw.gz, mw_dust.raw.gz and their metadata, and detail.raw.gz, a cube of noise
for clumps finer than the voxels; with a third folder, also pictures of
the layers read from the picture and each layer of light as a grid of its own, for tuning.
"""
import gzip
import sys

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates

PICTURE, OUT = sys.argv[1], sys.argv[2]
TUNING = sys.argv[3] if len(sys.argv) > 3 else None
Image.MAX_IMAGE_PIXELS = None

# Where the Sun and the galactic centre are in the picture, measured on its labelled version:
# the Sun at the centre of the distance rings, the centre at the brightest point of the bar.
# Longitude 0 is up and longitude 90 to the left.
SUN_PX = (2800.0, 3874.0)
CENTRE_PX = (2800.0, 2789.0)
SUN = 8.0                                           # kpc from the Sun to the centre, as in the map
PC_PER_PX = SUN * 1000.0 / (SUN_PX[1] - CENTRE_PX[1])

N_PLANE = 640                                       # cells along each axis of the plane
HALF = 20.0                                         # kpc
LIGHT_N, LIGHT_HALF = 96, 2.0                       # cells and kpc, above and below the plane
DUST_N, DUST_HALF = 32, 0.5
GAMMA = 2.0                                         # stored value = density ** (1 / GAMMA)
# The picture is drawn for the eye: next to its disc, its bulge is fainter than the real one is
# in photographs of the sky. This is how much brighter the bulge is made here.
BULGE_GAIN = 14.0
LOCAL_HOLE = (1.05, 1.25)
# The optical depth of the measured dust around the Sun, seen straight from above: 0.22 in a
# typical place (the median; the mean is 0.29, the spread between cells of this grid 0.8 of the
# mean). The spread-out dust of the model is as deep around there, so that the measured part
# does not stand out from afar. Measured dust is everywhere, not only in lanes.
MEASURED_COLUMN = 0.22                           # kpc from the Sun: no dust of this model inside
CELL = 2 * HALF / N_PLANE                           # kpc


def axis(n, half):
    return (np.arange(n, dtype=np.float32) + 0.5) / n * 2 * half - half


def smoothstep(a, b, v):
    t = np.clip((v - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


# ---- the picture, on the cells of the plane
# First axis: from the centre away from the Sun. Third axis: towards longitude 90.
x = axis(N_PLANE, HALF)[:, None]
z = axis(N_PLANE, HALF)[None, :]
col = SUN_PX[0] - z * 1000.0 / PC_PER_PX + 0 * x
row = SUN_PX[1] - (x + SUN) * 1000.0 / PC_PER_PX + 0 * z
picture = (np.asarray(Image.open(PICTURE).convert("RGB"), np.float32) / 255.0) ** 2.2
blur = CELL * 1000.0 / PC_PER_PX / 2.5              # so that a cell averages the pixels it covers
plane = np.stack([map_coordinates(gaussian_filter(picture[..., c], blur), [row, col], order=1, cval=0.0) for c in range(3)], axis=-1)
del picture
r = np.sqrt(x * x + z * z)
R, G, B = plane[..., 0], plane[..., 1], plane[..., 2]
lum = plane.mean(axis=-1)


def wide(v, kpc):
    return gaussian_filter(v, kpc / CELL)


# The bar and bulge: the yellow, bright middle.
bulge = smoothstep(0.8, 1.05, wide(R, 0.3) / (wide(B, 0.3) + 1e-4)) * smoothstep(0.05, 0.15, wide(lum, 0.3)) * smoothstep(6.0, 4.5, r)
bulge = wide(bulge, 0.45)
# Dust lanes: darker than their surroundings and not blue.
around = wide(lum, 0.25)
lanes = np.clip(1.0 - lum / (around + 1e-4), 0.0, 1.0) * smoothstep(0.7, 0.95, R / (B + 1e-4)) * (1.0 - bulge)
lanes = np.clip(lanes * 2.5, 0.0, 1.0)
# The light without the lanes drawn over it: under a lane, what is around it.
keep = 1.0 - smoothstep(0.1, 0.5, lanes)
filled = np.stack([gaussian_filter(plane[..., c] * keep, 3.0) / (gaussian_filter(keep, 3.0) + 1e-3) for c in range(3)], axis=-1)
shine = plane * keep[..., None] + filled * (1.0 - keep[..., None])
# Three layers of light: the smooth haze of old stars, what stands out of it (arms, clusters,
# gas clouds), and the bulge.
haze = np.minimum(shine, np.stack([wide(shine[..., c], 0.6) for c in range(3)], axis=-1))
arms = shine - haze
haze_disc, arms_disc, bulge_light = haze * (1 - bulge[..., None]), arms * (1 - bulge[..., None]), shine * bulge[..., None]


def layer(n, half, profile):
    """A vertical profile on the cells of a grid, with an integral of 1 over the height (in kpc)."""
    y = axis(n, half)
    p = profile(np.abs(y)) * smoothstep(half, half * 0.8, np.abs(y))
    return (p / (p.sum() * 2 * half / n)).astype(np.float32)


def write(name, density, rgb, top, n, half, note, folder=OUT):
    """density: (N, n, N) floats; rgb: (N, n, N, 3) or (3,) colors in 0..1."""
    out = np.empty(density.shape + (4,), np.uint8)
    out[..., :3] = (np.clip(rgb, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    out[..., 3] = (np.clip(density / top, 0.0, 1.0) ** (1.0 / GAMMA) * 255.0 + 0.5).astype(np.uint8)
    with gzip.open(f"{folder}/{name}.raw.gz", "wb", compresslevel=4) as f:
        f.write(out.tobytes())                                  # the third axis runs fastest
    with open(f"{folder}/{name}_metadata.txt", "w") as f:
        f.write("# Made by tools/make_galaxy.py from the face-on map of the Milky Way by R. Hurt (NASA/JPL-Caltech/ESO).\n")
        f.write("# Axes: (from the galactic centre away from the Sun, north, longitude 90).\n")
        f.write(f"Grid size: {N_PLANE}x{n}x{N_PLANE}\n")
        f.write(f"Volume bounds: [{-HALF}, {-half}, {-HALF}] to [{HALF}, {half}, {HALF}]\n")
        f.write(f"# Bounds are in kpc: sizePc 1000, centrePc [0, 0, {SUN * 1000:.0f}], gamma {GAMMA:.0f}. {note}\n")


# ---- light: each layer of the plane times how it is spread in height
thin = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / 0.09))
thick = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / 0.30))
round_ = layer(LIGHT_N, LIGHT_HALF, lambda y: np.exp(-y / 0.5))
# The bulge is as tall as it is deep inside its outline in the plane, so that it is round from the side too.
reach = wide(bulge, 0.7)
reach = np.maximum(1.9 * reach / reach.max(), 1e-3)
rounded = smoothstep(reach[:, None, :], 0.45 * reach[:, None, :], np.abs(axis(LIGHT_N, LIGHT_HALF))[None, :, None] + 0 * reach[:, None, :])
light = (haze_disc[:, None, :, :] * thick[None, :, None, None] + arms_disc[:, None, :, :] * thin[None, :, None, None]
         + BULGE_GAIN * bulge_light[:, None, :, :] * round_[None, :, None, None] * rounded[..., None])
if TUNING:                                           # each layer of light in a grid of its own
    for name, part in (("haze", haze_disc[:, None, :, :] * thick[None, :, None, None]), ("arms", arms_disc[:, None, :, :] * thin[None, :, None, None]),
                       ("bulge", bulge_light[:, None, :, :] * round_[None, :, None, None])):
        s = part.max(axis=-1)
        top = float(np.quantile(s[s > 0], 0.9995))
        write(f"layer_{name}", s, part / np.maximum(s, 1e-9)[..., None], top, LIGHT_N, LIGHT_HALF, f"top {top:.6g}", TUNING)
        print(f"layer {name}: top {top:.6g}")
strength = light.max(axis=-1)
light_top = float(np.quantile(strength[strength > 0], 0.9995))
write("mw_light", strength, light / np.maximum(strength, 1e-9)[..., None], light_top, LIGHT_N, LIGHT_HALF,
      f"Seen from above, lightPerKpc 1 gives the picture (its bulge {BULGE_GAIN:.0f} times brighter); full density is {light_top:.4g} of the picture's white per kpc.")
print(f"light: full density {light_top:.3g} per kpc")
del light, strength

# ---- dust: optical depth per kpc in visible light
rng = np.random.default_rng(7)


def fractal(shape, slope, stretch):
    """Noise with structure on all scales, mean 0 and deviation 1."""
    white = np.fft.rfftn(rng.standard_normal(shape, dtype=np.float32))
    k = np.sqrt(sum(np.meshgrid(np.fft.fftfreq(shape[0]) ** 2, (np.fft.fftfreq(shape[1]) * stretch) ** 2, np.fft.rfftfreq(shape[2]) ** 2, indexing="ij")))
    k[0, 0, 0] = 1.0
    out = np.fft.irfftn(white * k ** (-slope), s=shape).astype(np.float32)
    return (out - out.mean()) / out.std()


shape = (N_PLANE, DUST_N, N_PLANE)
# Everywhere in the disc: about a magnitude per kpc at the Sun, more inwards and in the arms.
spread = 1.1 * np.exp(-(np.clip(r, 3.0, None) - SUN) / 4.5) * smoothstep(18.0, 15.0, r) * (0.4 + 0.6 * smoothstep(1.5, 3.5, r))
spread = spread * (1.0 + 1.5 * np.clip(wide(arms_disc.mean(axis=-1), 0.2) / 0.05, 0.0, 1.0))
height = layer(DUST_N, DUST_HALF, lambda y: np.exp(-y / 0.09)) * 0.18           # 1 in the plane
y3 = axis(DUST_N, DUST_HALF)[None, :, None]
from_sun = np.sqrt((x[:, None, :] + SUN) ** 2 + y3 ** 2 + z[:, None, :] ** 2)
# As much dust around the measured part as in it: the level of the spread-out dust is set so
# that a typical place there is as deep as a typical place of the measured map.
around = (np.sqrt((x + SUN) ** 2 + z ** 2) > 1.4) & (np.sqrt((x + SUN) ** 2 + z ** 2) < 3.0)
# Clumps as uneven as the measured ones.
clumps = np.exp(0.8 * fractal(shape, 1.4, N_PLANE / DUST_N * DUST_HALF / HALF * 0.5) - 0.32)
level = 1.0
for _ in range(3):
    dust = (level * spread + 10.0 * lanes)[:, None, :] * height[None, :, None] * clumps
    level *= MEASURED_COLUMN / np.median((dust.sum(axis=1) * 2 * DUST_HALF / DUST_N)[around])
print(f"dust: the spread-out part is {level:.2f} times its first guess, to match the measured column")
dust *= smoothstep(LOCAL_HOLE[0], LOCAL_HOLE[1], from_sun)
dust_top = float(np.quantile(dust[dust > 0], 0.9995))
write("mw_dust", dust, np.array([0.55, 0.38, 0.24], np.float32), dust_top, DUST_N, DUST_HALF, f"dustPerKpc {dust_top:.1f}.")
through = dust.sum(axis=1) * 2 * DUST_HALF / DUST_N
print(f"dust: full density {dust_top:.3g} per kpc; seen from above the lanes have a depth of {np.quantile(through[lanes > 0.5], 0.5):.2f}, the rest {np.median(through[(r > 4) & (r < 12) & (lanes < 0.05)]):.2f}; "
      f"around the measured part: typical {np.median(through[around]):.2f}, mean {through[around].mean():.2f}, spread {through[around].std() / through[around].mean():.2f} of it (measured: {MEASURED_COLUMN}, 0.29, 0.8)")
print(f'for the map: light "lightPerKpc" in units of {light_top:.4g}; dust "dustPerKpc": {dust_top:.1f}')


# ---- noise for the clumps finer than the voxels: a cube that repeats, a bell-shaped spread
# around the byte 128 with 48 to one deviation, with structure on all its scales.
DETAIL = 64
g = np.random.default_rng(23)
white = np.fft.rfftn(g.standard_normal((DETAIL,) * 3, dtype=np.float32))
k = np.sqrt(sum(np.meshgrid(np.fft.fftfreq(DETAIL) ** 2, np.fft.fftfreq(DETAIL) ** 2, np.fft.rfftfreq(DETAIL) ** 2, indexing="ij")))
k[0, 0, 0] = 1.0
field = np.fft.irfftn(white * k ** (-1.3), s=(DETAIL,) * 3)
field = (field - field.mean()) / field.std()
with gzip.open(f"{OUT}/detail.raw.gz", "wb") as f:
    f.write(np.clip(128.0 + 48.0 * field + 0.5, 0, 255).astype(np.uint8).tobytes())

if TUNING:
    for name, v in (("light", shine), ("bulge", bulge), ("lanes", lanes), ("dust_from_above", 1.0 - np.exp(-through))):
        Image.fromarray((np.clip(v[::-1], 0, 1) ** (1 / 2.2) * 255).astype(np.uint8)).save(f"{TUNING}/preview_{name}.png")
