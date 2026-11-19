#!/usr/bin/env python3
"""A galaxy's measured light along its long axis, for tools/galaxy_xyz.py --profile: its surface
brightness (mag/arcsec²) at each distance from its middle (kpc), from a calibrated survey image.

The pictures the image model makes are made to be looked at: a galaxy's middle, hundreds of times
brighter than its disc, is as white as the arms. galaxy_xyz.py puts the light back from a table of
the real one; this measures that table for any galaxy in the DESI Legacy Surveys' sky (about
two thirds of it, dec -70 to +85, not near the Milky Way's plane), from its r-band image, whose
pixels are calibrated: nanomaggies, a magnitude of 22.5 - 2.5 log10(flux).

  middle     the peak of the image, smoothed, near the asked place;
  long axis  from the second moments of the brightest 15% of it (its disc, not only its bulge);
  rings      ellipses of the given axis ratio (b/a, as galaxy_fit.py measures it), each the
             median of its pixels (stars in it do not count), less the sky (the median of the
             outermost ones), out to where it is fainter than 26.5 mag/arcsec².
The innermost arcsecond is smoothed by the seeing (about 1.2"): the true middle is brighter.

Usage: galaxy_profile.py OUT.json --ra DEG --dec DEG --distance-pc D --axis-ratio Q --tilt DEG
                         [--size-arcmin 16] [--band r] [--name NAME]
Writes {"tiltDeg", "source", "table": [[kpc, mag/arcsec²], ...]} (galaxy_xyz.py --profile OUT.json).
"""
import argparse
import json
import urllib.request

import numpy as np
from scipy.ndimage import gaussian_filter

args = argparse.ArgumentParser()
args.add_argument("out")
args.add_argument("--ra", type=float, required=True)
args.add_argument("--dec", type=float, required=True)
args.add_argument("--distance-pc", type=float, required=True)
args.add_argument("--axis-ratio", type=float, required=True, help="b/a of its disc as seen (galaxy_fit.py's photoMetadata.axisRatio)")
args.add_argument("--tilt", type=float, required=True, help="its tilt, degrees (fitted.json tiltDeg): galaxy_xyz.py makes the disc face on with it")
args.add_argument("--size-arcmin", type=float, default=16.0)
args.add_argument("--band", default="r")
args.add_argument("--name", default="")
args.add_argument("--fits", help="a FITS file already fetched (else from legacysurvey.org)")
a = args.parse_args()

PS = 0.5  # arcsec a pixel asked for


def read_fits(raw):
    head, i = {}, 0
    while True:
        card = raw[i:i + 80].decode("ascii", "replace")
        i += 80
        if card[:8].strip() == "END":
            break
        if card[8:10] == "= ":
            head[card[:8].strip()] = card[10:].split("/")[0].strip().strip("'").strip()
    i = ((i + 2879) // 2880) * 2880
    nx, ny = int(head["NAXIS1"]), int(head["NAXIS2"])
    return np.frombuffer(raw[i:i + nx * ny * 4], ">f4").reshape(ny, nx).astype(np.float64)


size = int(a.size_arcmin * 60 / PS)
if a.fits:
    raw = open(a.fits, "rb").read()
    source = a.fits
else:
    url = (f"https://www.legacysurvey.org/viewer/cutout.fits?ra={a.ra}&dec={a.dec}&layer=ls-dr10"
           f"&pixscale={PS}&size={size}&bands={a.band}")
    raw = urllib.request.urlopen(url, timeout=600).read()
    source = url
img = read_fits(raw)
bad = ~np.isfinite(img)
img[bad] = 0.0
H, W = img.shape

smooth = gaussian_filter(img, 4)
box = 100
cy, cx = np.unravel_index(smooth[H // 2 - box:H // 2 + box, W // 2 - box:W // 2 + box].argmax(), (2 * box, 2 * box))
cy, cx = cy + H // 2 - box, cx + W // 2 - box
ys, xs = np.nonzero(smooth > np.quantile(smooth, 0.85))
values, vectors = np.linalg.eigh(np.cov(np.stack([xs - cx, ys - cy])))
L = vectors[:, 1]
S = np.array([-L[1], L[0]])
yy, xx = np.mgrid[0:H, 0:W]
u = (xx - cx) * L[0] + (yy - cy) * L[1]
v = ((xx - cx) * S[0] + (yy - cy) * S[1]) / a.axis_ratio
r = np.hypot(u, v) * PS
pc_per_arcsec = a.distance_pc * np.pi / 180 / 3600
sky = float(np.median(img[(r > 0.9 * r.max()) & ~bad]))

edges = np.unique(np.concatenate([np.arange(0, 10, 1.0), np.geomspace(10, r.max(), 40)]))
table = []
for lo, hi in zip(edges[:-1], edges[1:]):
    k = (r >= lo) & (r < hi) & ~bad
    if k.sum() < 5:
        continue
    f = float(np.median(img[k])) - sky
    if f <= 0:
        break
    mu = 22.5 - 2.5 * np.log10(f / PS ** 2)
    table.append([round((lo + hi) / 2 * pc_per_arcsec / 1000, 4), round(mu, 2)])
    if mu > 26.5:
        break

out = {"name": a.name, "tiltDeg": a.tilt,
       "source": f"DESI Legacy Surveys DR10 {a.band} band ({source}), median in ellipses of b/a {a.axis_ratio} along the long axis "
                 f"(tools/galaxy_profile.py); the innermost arcsecond smoothed by the seeing",
       "table": table}
json.dump(out, open(a.out, "w"), indent=1)
print(f"{a.out}: {len(table)} rings, {table[0][1]:.2f} mag/arcsec² in the middle to {table[-1][1]:.2f} at {table[-1][0]:.1f} kpc "
      f"(sky left {sky:.4f} nanomaggies a pixel; the bright body's b/a {np.sqrt(values[0] / values[1]):.2f})")
