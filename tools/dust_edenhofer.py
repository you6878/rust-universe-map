"""Turns the measured dust map of Edenhofer et al. (2023) into a voxel grid for the viewer.

The map gives the dust density around the Sun, from 69 pc to 1.25 kpc, on HEALPix spheres at
516 distances (https://doi.org/10.5281/zenodo.8187943, CC BY 4.0). fetch_shells.sh downloads
the spheres one file each; any subset of them works here, the more the finer in distance.

Usage: dust_edenhofer.py SHELL_DIR RADII.npy OUT_DIR
Writes local_dust_measured.raw.gz and its metadata. match_panorama.py makes the grid of the map from it.
"""
import glob
import gzip
import re
import sys

import healpy as hp
import numpy as np

SHELLS, RADII, OUT = sys.argv[1:4]
N = (400, 128, 400)                 # voxels along (galactic centre, north, longitude 90)
HALF = (1250.0, 400.0, 1250.0)      # parsecs
FADE = (1050.0, 1245.0)             # the map hands over to the model of the galaxy between these
GAMMA = 2.0
# One unit of the map's extinction E is 2.8 magnitudes in V, and a magnitude is 0.921 of optical depth.
DEPTH_PER_E = 2.8 * 0.921
NPIX = 12 * 256 * 256

radii = np.load(RADII)
files = sorted(f for f in glob.glob(f"{SHELLS}/shell_*.bin") if hp.nside2npix(256) * 4 == __import__("os").path.getsize(f))
index = np.array([int(re.search(r"shell_(\d+)", f).group(1)) for f in files])
if len(files) < 2:
    sys.exit("need at least two complete shells")
shells = np.stack([np.fromfile(f, ">f4").astype(np.float32) for f in files])
log_r = np.log(radii[index])
print(f"{len(files)} shells, {radii[index[0]]:.0f} to {radii[index[-1]]:.0f} pc")


def axis(n, half):
    return (np.arange(n, dtype=np.float32) + 0.5) / n * 2 * half - half


y = axis(N[1], HALF[1])[None, :, None]
z = axis(N[2], HALF[2])[None, None, :]
grid = np.zeros(N, np.float32)
for i, x in enumerate(axis(N[0], HALF[0])):          # one slab at a time, to keep the memory small
    r = np.sqrt(x * x + y * y + z * z)
    lon = np.arctan2(z, x) + 0.0 * r                  # longitude grows from the centre towards l = 90
    colat = np.arccos(np.clip(y / np.maximum(r, 1e-6), -1.0, 1.0))
    pix = hp.ang2pix(256, colat, lon, nest=True)
    # Between which two shells the voxel lies, and how far along in log distance.
    lr = np.log(np.maximum(r, 1e-6))
    hi = np.clip(np.searchsorted(log_r, lr), 1, len(log_r) - 1)
    t = np.clip((lr - log_r[hi - 1]) / (log_r[hi] - log_r[hi - 1]), 0.0, 1.0)
    value = shells[hi - 1, pix] * (1.0 - t) + shells[hi, pix] * t
    value[(r < radii[index[0]]) | (r > radii[index[-1]])] = 0.0
    grid[i] = value * (1.0 - np.clip((r - FADE[0]) / (FADE[1] - FADE[0]), 0.0, 1.0))

top = float(np.quantile(grid[grid > 0], 0.9995))
out = np.empty(N + (4,), np.uint8)
out[..., :3] = np.array([140, 97, 61], np.uint8)     # the color the dust glows in
out[..., 3] = (np.clip(grid / top, 0.0, 1.0) ** (1.0 / GAMMA) * 255.0 + 0.5).astype(np.uint8)
with gzip.open(f"{OUT}/local_dust_measured.raw.gz", "wb", compresslevel=4) as f:
    f.write(out.tobytes())
with open(f"{OUT}/local_dust_measured_metadata.txt", "w") as f:
    f.write("# Dust around the Sun from Edenhofer et al. (2023), https://doi.org/10.5281/zenodo.8187943, CC BY 4.0.\n")
    f.write(f"# Made by tools/dust_edenhofer.py from {len(files)} of the 516 distance shells. Axes: (galactic centre, north, longitude 90).\n")
    f.write(f"Grid size: {N[0]}x{N[1]}x{N[2]}\n")
    f.write(f"Volume bounds: [{-HALF[0]}, {-HALF[1]}, {-HALF[2]}] to [{HALF[0]}, {HALF[1]}, {HALF[2]}]\n")
    f.write(f"# Bounds are in parsecs: sizePc 1, centrePc [0, 0, 0], gamma {GAMMA:.0f}, dustPerKpc {DEPTH_PER_E * top * 1000:.1f}.\n")
plane = grid[:, N[1] // 2 - 4:N[1] // 2 + 4, :]
print(f"full density {top:.4g} E/pc; in the plane the mean is {plane.mean():.4g} E/pc = {plane.mean() * 2.8 * 1000:.2f} mag/kpc in V")
print(f'for the map: "sizePc": 1, "centrePc": [0, 0, 0], "gamma": {GAMMA:.0f}, "dustPerKpc": {DEPTH_PER_E * top * 1000:.1f}')
