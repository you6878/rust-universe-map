#!/usr/bin/env python3
"""A galaxy's photograph as its disc from afar is drawn with (godot/scripts/galaxy_look.gd): the
galaxy alone, the foreground stars of our own galaxy and the sky around taken out.

Drawn on the galaxy's disc, a star of the photograph would be as large as a star cluster of the
galaxy, a smear thousands of light years across. So:
  body   the galaxy: where the photograph, much blurred, is brighter than the sky by a share of its
         middle, the piece that holds the middle, grown a little and with a soft edge; outside it
         black;
  stars  in the body, what stands out of the light around it (a --median pixel median) by far
         and is small, and --grow pixels round it: the median there instead. A smooth galaxy (M65)
         takes a wide median; a knotty one (M101: its HII regions are its own) a narrow one.
The photograph stays as it is for fitting (tools/galaxy_fit.py, galaxy_xyz.py): this is only what
is shown.

Usage: photo_clean.py PHOTO OUT [--body 0.02] [--stars 6] [--median 13] [--grow 3]
"""
import argparse

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, label, median_filter, binary_dilation

args = argparse.ArgumentParser()
args.add_argument("photo")
args.add_argument("out")
args.add_argument("--body", type=float, default=0.02, help="the body: brighter than the sky by this share of its middle")
args.add_argument("--stars", type=float, default=6.0, help="a star: above the median around it by this many times the noise")
args.add_argument("--median", type=int, default=13, help="the light around a star: a median this many pixels across (of 1024)")
args.add_argument("--grow", type=int, default=3, help="and this many pixels round each star taken out too")
a = args.parse_args()

rgb = np.asarray(Image.open(a.photo).convert("RGB"), np.float64) / 255.0
H, W = rgb.shape[:2]
lum = rgb.mean(axis=-1)
sky = np.median(lum)
noise = 1.4826 * np.median(np.abs(lum - sky))

# the body
smooth = gaussian_filter(lum, 6.0 * W / 1024)
cy, cx = np.unravel_index(gaussian_filter(lum, 20.0 * W / 1024).argmax(), lum.shape)
pieces, _ = label(smooth > sky + a.body * (smooth[cy, cx] - sky))
body = pieces == pieces[cy, cx]
body = binary_dilation(body, iterations=int(12 * W / 1024))
soft = np.clip(gaussian_filter(body.astype(np.float64), 8.0 * W / 1024) * 2.0 - 0.5, 0.0, 1.0)

# the stars in it
k = int(a.median * W / 1024) | 1
around = np.stack([median_filter(rgb[..., c], size=k) for c in range(3)], -1)
excess = lum - around.mean(axis=-1)
star = excess > a.stars * max(noise, 1e-3)
# small: not the middle of the galaxy itself
star &= gaussian_filter(lum, 20.0 * W / 1024) < 0.9 * gaussian_filter(lum, 20.0 * W / 1024).max()
star = binary_dilation(star, iterations=max(1, int(a.grow * W / 1024)))
clean = np.where(star[..., None], around, rgb)
clean = np.clip((clean - sky) / max(1.0 - sky, 1e-6), 0.0, 1.0) * soft[..., None]

Image.fromarray((clean * 255.0 + 0.5).astype(np.uint8)).save(a.out, quality=95)
print(f"{a.out}: the body {body.mean() * 100:.1f}% of the picture, {int(label(star)[1])} stars taken out")
