"""Puts a galaxy's three views (one along each axis, made by the image model) on one grid: the same
size and the same parsecs a pixel in all three, and finds where that grid lies in the galaxy.

The views' own axes: U down the face view's rows, W along its columns (both in the disc's plane),
Y its pole (up).
  y  seen along Y, from above: the face view as it is; rows +U, columns +W. It sets the scale.
  x  seen along U (the image model drew it from the face view's bottom edge): columns +W, rows -Y.
  z  seen along W (drawn from the face view's right edge, its top on the left): columns +U, rows -Y.
The face view was made from the photograph stretched face on, in the photograph's own orientation:
its rows and columns are not the disc's long and short axes. With --photo (fitted.json and the
galaxy's name) the face view is laid on the disc as the photograph sees it (fitted.json's
photoMetadata: a place X along the long axis, Z along the short one, is at W = X lx + Z sx,
U = X ly + Z sy, (lx, ly) and (sx, sy) the axes' directions in the photograph), and the parsecs a
pixel and the disc's middle on the face view are searched for: those that make it, seen from the
Earth, most like the photograph.

The image model draws the edge-on views at a scale of its own. Each is scaled about its middle so
that the width holding 99% of its light is the face view's along the same axis, and moved so that
its middle (the peak of its light, the bulge) is the face view's along that axis, and the mid-plane is the
middle row. Prints the widths before and after.

The edge views' mid-plane is put on the middle row: where their dust lane is (the mid-plane of a
disc), else the peak of their light.

The pictures are made to be looked at: their light is squeezed (Andromeda's middle, hundreds of times
brighter than its disc, is as white as the arms 3 kpc out). With --profile, a table of the galaxy's
measured surface brightness along its long axis (kpc, mag/arcsec^2), the light is put back: a curve
from a pixel's light to its true light, found where the face view is not white (the mean of a ring
of the face view against the table's brightness at that radius, face on); where it is white the
table's brightness at its radius; in the edge views, where they are white, the light below and above
followed on into it (log-linearly) to a peak. Written as 16-bit pictures of log10 of the light
(OUT_PREFIX_{y,x,z}_light.png), the 8-bit ones kept for their colours.

Usage: galaxy_xyz.py FACE ALONG_U ALONG_W OUT_PREFIX [--size 1024] [--stars 1e12]
                      [--photo FITTED.json PHOTO.jpg] [--name NAME] [--profile M31|PROFILE.json]
Writes OUT_PREFIX_y.png, _x.png, _z.png (size x size RGB) and OUT_PREFIX.json, the frame.
"""
import argparse

import numpy as np
from scipy.ndimage import gaussian_filter1d
from PIL import Image

args = argparse.ArgumentParser()
args.add_argument("face")
args.add_argument("along_u")
args.add_argument("along_w")
args.add_argument("out")
args.add_argument("--size", type=int, default=1024)
args.add_argument("--half-kpc", type=float, default=39.863)
args.add_argument("--stars", type=float, default=1e12)
args.add_argument("--name", default="Andromeda Galaxy")
args.add_argument("--photo", nargs=2, metavar=("FITTED", "PHOTO"))
args.add_argument("--profile", help="M31 (the table below), or a JSON file of tools/galaxy_profile.py: {tiltDeg, source, table}")
args.add_argument("--middle", nargs=2, type=float, metavar=("COL", "ROW"), help="the galaxy's middle on the face view, where it is known (the Milky Way's map is cut round its centre): not its brightest place")
args.add_argument("--height-factor", type=float, default=1.0, help="the edge views made this many times taller about their mid-plane, after they are laid on the face view: their thickness set to the measured one (the Milky Way's)")
args.add_argument("--kpc-range", nargs=2, type=float, metavar=("LO", "HI"), help="the kpc a pixel searched for with --photo (default 0.01 0.2, Andromeda's)")

# Surface brightness along the long axis, R band, mag/arcsec^2 at kpc from the middle: Andromeda's,
# read off Courteau et al. 2011 (ApJ 739, 20, figure of the R-band major-axis profile), to about
# half a magnitude. Seen tilted (77 deg): face on, the disc is fainter by 2.5 log10(1 / cos i).
PROFILES = {"M31": {"tiltDeg": 77.0, "source": "Courteau et al. 2011, ApJ 739, 20: R-band major-axis profile, read to ~0.5 mag",
                    "table": [[0.05, 14.5], [0.1, 15.2], [0.25, 16.0], [0.5, 16.8], [1, 17.6], [2, 18.6], [3, 19.1], [5, 19.8],
                              [8, 20.4], [10, 20.8], [12, 21.2], [15, 21.7], [20, 22.7], [25, 23.8], [30, 25.0], [35, 26.2]]}}
a = args.parse_args()
N = a.size


def load(f):
    return np.asarray(Image.open(f).convert("RGB").resize((N, N), Image.LANCZOS), dtype=np.float64) / 255.0


def light(img):
    return (img ** 2.2).mean(axis=2)


def extent(profile):
    """Where the 0.5% and 99.5% points of a light profile are, and its middle: the peak of the
    profile smoothed over 2% of the picture (the bulge; the centroid of all the light can be kpc
    away from it in a lopsided disc) (pixels)."""
    p = profile / profile.sum()
    c = np.cumsum(p)
    smooth = gaussian_filter1d(p, N * 0.02)
    k = int(np.argmax(smooth))
    # the peak to a fraction of a pixel: a parabola through it and its neighbours
    if 0 < k < len(p) - 1:
        a, b, d = smooth[k - 1], smooth[k], smooth[k + 1]
        k = k + 0.5 * (a - d) / (a - 2 * b + d) if a - 2 * b + d < 0 else k
    return np.searchsorted(c, 0.005), np.searchsorted(c, 0.995), float(k)


def dust_lane(l):
    """The row of an edge-on view's dust lane: in each column of the disc's middle half, the darkest
    row between two brighter ones near the middle row; their median. None if there is none."""
    rows = []
    cols = np.nonzero(l.sum(axis=0) > 0.05 * l.sum(axis=0).max())[0]
    if len(cols) == 0:
        return None
    lo, hi = np.quantile(cols, 0.2), np.quantile(cols, 0.8)
    # within the disc's light: around the brightest row of the middle half
    m = int(np.argmax(l[:, int(lo):int(hi)].sum(axis=1)))
    w = 12
    for c in range(int(lo), int(hi), 2):
        p = l[m - w:m + w + 1, c]
        k = int(np.argmin(p[2:-2])) + 2
        if p[k] < 0.5 * min(p[:k].max(), p[k + 1:].max()):
            rows.append(m - w + k)
    return float(np.median(rows)) if len(rows) > 10 else None


def place(img, scale, dx, dy):
    """The picture scaled by `scale` about its middle, then moved by (dx, dy) pixels."""
    h, w = img.shape[:2]
    m = (N - 1) / 2.0
    # each output pixel reads the input at (out - m - d) / scale + m
    out = Image.fromarray((np.clip(img, 0, 1) * 255 + 0.5).astype(np.uint8))
    inv = 1.0 / scale
    coeffs = (inv, 0, m - (m + dx) * inv, 0, inv, m - (m + dy) * inv)
    return np.asarray(out.transform((N, N), Image.AFFINE, coeffs, resample=Image.BICUBIC), dtype=np.float64) / 255.0


face = load(a.face)
f = light(face)
fx = extent(f.sum(axis=0))      # along the face view's horizontal axis
fy = extent(f.sum(axis=1))      # ... and its vertical one
if a.middle:
    fx = (fx[0], fx[1], a.middle[0])
    fy = (fy[0], fy[1], a.middle[1])
results = {"y": face}
for name, file, target in (("x", a.along_u, fx), ("z", a.along_w, fy)):
    img = load(file)
    l = light(img)
    ex = extent(l.sum(axis=0))
    ey = extent(l.sum(axis=1))
    scale = (target[1] - target[0]) / (ex[1] - ex[0])
    m = (N - 1) / 2.0
    # after scaling about the middle, the centroid is at m + (c - m) * scale: move it to the target
    dx = target[2] - (m + (ex[2] - m) * scale)
    dy = m - (m + (ey[2] - m) * scale)
    placed = place(img, scale, dx, dy)
    # the mid-plane: the dust lane, where it crosses the disc's middle half
    lane = dust_lane(light(placed))
    if lane is not None:
        placed = place(placed, 1.0, 0.0, m - lane)
    if a.height_factor != 1.0:
        # taller (or flatter) about the mid-plane, the middle row
        im = Image.fromarray((np.clip(placed, 0, 1) * 255 + 0.5).astype(np.uint8))
        k = 1.0 / a.height_factor
        placed = np.asarray(im.transform((N, N), Image.AFFINE, (1, 0, 0, 0, k, m - m * k), resample=Image.BICUBIC), dtype=np.float64) / 255.0
    pl = light(placed)
    after = extent(pl.sum(axis=0))
    print(f"{name}: width {ex[1] - ex[0]} px -> {after[1] - after[0]} px (face {target[1] - target[0]} px), scale {scale:.4f}, "
          f"middle {ex[2]:.1f} -> {after[2]:.1f} (face {target[2]:.1f}), dust lane at row {lane if lane is None else round(lane, 1)} -> {(N - 1) / 2:.1f}")
    results[name] = placed
import json
import os

# ---- where the grid lies in the galaxy: as the photograph sees the face view
frame = {"kpcPerPx": 2 * a.half_kpc / N, "middlePx": {"col": (N - 1) / 2, "row": (N - 1) / 2}, "wAxis": [1.0, 0.0], "uAxis": [0.0, 1.0], "likePhotograph": None}
if a.photo:
    from scipy.ndimage import map_coordinates
    d = json.load(open(a.photo[0]))
    entry = next(g for g in (d if isinstance(d, list) else d.get("galaxies", d)) if g["name"] == a.name)
    m = entry["photoMetadata"]
    photo = (np.asarray(Image.open(a.photo[1]).convert("RGB"), np.float32) / 255.0) ** 2.2
    photo = photo.mean(-1)
    edge = np.concatenate([photo[0], photo[-1], photo[:, 0], photo[:, -1]])
    photo = np.clip(photo - np.median(edge) * 1.5, 0.0, None)
    H, W = photo.shape
    (pcx, pcy), kpp = m["middlePx"], m["kpcPerPx"]
    (lx, ly), (sx, sy), ci = m["longAxisPx"], m["shortAxisPx"], m["cosTilt"]
    rr, cc = np.mgrid[0:H:4, 0:W:4]
    dx, dy = cc - pcx, rr - pcy
    X = (dx * lx + dy * ly) * kpp           # where on the disc each pixel of the photograph is, kpc
    Z = (dx * sx + dy * sy) * kpp / ci
    target = photo[::4, ::4]
    top = light(face)

    def like(kf, c0, r0):
        s = map_coordinates(top, [r0 + (X * ly + Z * sy) / kf, c0 + (X * lx + Z * sx) / kf], order=1, cval=0.0)
        mm = (target > 0.02 * target.max()) | (s > 0.02 * s.max())
        return float(np.corrcoef(np.sqrt(s[mm]), np.sqrt(target[mm]))[0, 1])

    best = (-1.0,)
    if a.kpc_range:
        kfs, fine = np.geomspace(a.kpc_range[0], a.kpc_range[1], 80) * 1024 / N, 0.0
    else:
        kfs, fine = np.arange(0.01, 0.2, 0.0025) * 1024 / N, 0.002
    for kf in kfs:
        for c0 in np.arange(N * 0.47, N * 0.53, N / 256):
            for r0 in np.arange(N * 0.47, N * 0.53, N / 256):
                v = like(kf, c0, r0)
                if v > best[0]:
                    best = (v, kf, c0, r0)
    _, kf, c0, r0 = best
    step = fine or kf * 0.03
    for k2 in np.arange(kf - step, kf + step * 1.001, step / 10):
        for c in np.arange(c0 - N / 256, c0 + N / 256 + 0.1, 1):
            for r in np.arange(r0 - N / 256, r0 + N / 256 + 0.1, 1):
                v = like(k2, c, r)
                if v > best[0]:
                    best = (v, k2, c, r)
    v, kf, c0, r0 = best
    frame = {"kpcPerPx": float(kf), "middlePx": {"col": float(c0), "row": float(r0)}, "wAxis": [lx, sx], "uAxis": [ly, sy], "likePhotograph": v}
    print(f"seen from the Earth like the photograph {v:.3f}: {kf * 1000:.1f} pc a pixel, the disc's middle at column {c0:.1f}, row {r0:.1f}")
for name, img in results.items():
    Image.fromarray((np.clip(img, 0, 1) * 255 + 0.5).astype(np.uint8)).save(f"{a.out}_{name}.png")

# ---- the light put back (--profile)
light_info = None
if a.profile:
    prof = PROFILES[a.profile] if a.profile in PROFILES else json.load(open(a.profile))
    tab = np.array(prof["table"], dtype=np.float64)
    cos_i = np.cos(np.radians(prof["tiltDeg"]))
    # face on: the disc (beyond 3 kpc) fainter by 2.5 log10(1/cos i), the round bulge (inside 1 kpc) not
    ramp = np.clip((tab[:, 0] - 1.0) / 2.0, 0.0, 1.0)
    true_of_r = 10 ** (-0.4 * (tab[:, 1] + ramp * 2.5 * np.log10(1.0 / cos_i)))
    log_true = lambda r: np.interp(np.log10(np.maximum(r, tab[0, 0])), np.log10(tab[:, 0]), np.log10(true_of_r), right=None)
    pc = frame["kpcPerPx"]
    top = light(results["y"])
    rr = np.hypot(*np.meshgrid(np.arange(N) - frame["middlePx"]["col"], np.arange(N) - frame["middlePx"]["row"])) * pc
    WHITE = 0.9 * top.max()
    # the curve: a ring's mean light against the table's at its radius, where the ring is not white
    pts = []
    for r in np.arange(1.0, tab[-1, 0], 0.5):
        ring = (rr >= r - 0.25) & (rr < r + 0.25)
        mean = top[ring].mean()
        if 0 < mean < WHITE:
            pts.append((np.log10(mean), float(log_true(r))))
    pts = np.array(sorted(pts))
    # monotone: as the pixel's light grows its true light does not fall
    pts[:, 1] = np.maximum.accumulate(pts[:, 1])
    slope_lo = (pts[1, 1] - pts[0, 1]) / max(pts[1, 0] - pts[0, 0], 1e-9)

    def curve(l):
        x = np.log10(np.maximum(l, 1e-6))
        y = np.interp(x, pts[:, 0], pts[:, 1])
        return np.where(x < pts[0, 0], pts[0, 1] + (x - pts[0, 0]) * slope_lo, y)

    true_top = np.where(top < WHITE, curve(top), log_true(rr))
    true_top = np.where(top > 1e-5, true_top, -np.inf)
    if prof.get("byRadius"):
        # a picture whose light does not grow toward the middle as the real light does (the Milky
        # Way's map, drawn for the eye: its bar no brighter than its arms): each ring of it made
        # as bright, on average, as the table says, its arms and lanes kept within the ring
        bins = np.unique(np.concatenate([[0.0], np.geomspace(max(pc, 0.02), 1.0, 12), np.arange(1.0, rr.max() + 0.5, 0.5)]))
        idx = np.digitize(rr, bins)
        gain = np.zeros(len(bins) + 1)
        for k in range(1, len(bins)):
            ring = (idx == k) & (top > 1e-5)
            if ring.any():
                r_mid = 0.5 * (bins[k - 1] + bins[k])
                gain[k] = 10 ** float(log_true(r_mid)) / max(top[ring].mean(), 1e-9)
        true_top = np.where(top > 1e-5, np.log10(np.maximum(top * gain[idx], 1e-30)), -np.inf)
        curve = lambda l: np.log10(np.maximum(l, 1e-6))      # the edge views' heights as drawn

    def edge_true(img):
        l = light(img)
        t = np.where(l > 1e-5, curve(l), -np.inf)
        # where white: the light below and above followed on into it, log-linearly, to a peak
        for c in range(N):
            col = l[:, c]
            white = col >= WHITE
            if not white.any():
                continue
            k = 0
            while k < N:
                if not white[k]:
                    k += 1
                    continue
                e = k
                while e < N and white[e]:
                    e += 1
                lo_k, hi_k = k - 1, e        # the last ones not white, either side
                def side(at, step):
                    j = at - 2 * step
                    if 0 <= at < N and 0 <= j < N and np.isfinite(t[at, c]) and np.isfinite(t[j, c]):
                        return t[at, c], (t[at, c] - t[j, c]) / 2.0
                    return (t[at, c] if 0 <= at < N else -np.inf), 0.0
                v1, s1 = side(lo_k, -1)   # coming down from above: rising toward the white
                v2, s2 = side(hi_k, 1)
                for j in range(k, e):
                    a1 = v1 + max(s1, 0.0) * (j - lo_k)
                    a2 = v2 + max(s2, 0.0) * (hi_k - j)
                    t[j, c] = min(a1, a2) if np.isfinite(a1) and np.isfinite(a2) else max(a1, a2)
                k = e
        return t

    maps = {"y": true_top, "x": edge_true(results["x"]), "z": edge_true(results["z"])}
    finite = np.concatenate([v[np.isfinite(v)] for v in maps.values()])
    lo, hi = float(finite.min()), float(finite.max())
    for name, v in maps.items():
        q = np.where(np.isfinite(v), 1 + (v - lo) / (hi - lo) * 65534, 0)
        Image.fromarray(np.round(q).astype(np.uint16)).save(f"{a.out}_{name}_light.png")
    light_info = {"files": {k: f"{os.path.basename(a.out)}_{k}_light.png" for k in maps}, "log10Range": [lo, hi], "profile": prof,
                  "about": "the light put back (tools/galaxy_xyz.py --profile): 16-bit, 0 none, else log10(light) = lo + (v - 1) / 65534 (hi - lo)"}
    mid = (rr < 0.06)
    at8 = (rr > 7.8) & (rr < 8.2)
    print(f"light put back: {len(pts)} points of the curve; the middle {10 ** true_top[mid].mean() / np.mean(10 ** true_top[at8]):.0f} times as bright as 8 kpc out "
          f"(the picture: {top[mid].mean() / top[at8].mean():.1f} times)")
base = os.path.basename(a.out)
json.dump({
    "name": a.name, "size": N, "stars": a.stars,
    "views": {"y": f"{base}_y.png", "x": f"{base}_x.png", "z": f"{base}_z.png"},
    "frame": frame,
    "about": "a galaxy from three views (tools/galaxy_xyz.py), size x size pixels each, kpcPerPx a pixel: y seen from above, rows +U, "
             "columns +W; x seen along U, columns +W, rows -Y; z seen along W, columns +U, rows -Y (Y up, out of the disc). In the "
             "galaxy's own frame (fitted.json: X its long axis, Z its short one, parsecs from its middle) a place is at "
             "W = X wAxis[0] + Z wAxis[1], U = X uAxis[0] + Z uAxis[1], so at column middlePx.col + W / kpcPerPx, row middlePx.row + "
             "U / kpcPerPx of the views. Light is sRGB^2.2.",
    "bulgePx": {"col": fx[2], "row": fy[2]},
    **({"light": light_info} if light_info else {}),
}, open(f"{a.out}.json", "w"), indent=1)
print(f"middle (the bulge) in the face view: col {fx[2]:.1f}, row {fy[2]:.1f}; the mid-plane row {(N - 1) / 2:.1f}")
print(f"{N} x {N} each, {frame['kpcPerPx'] * 1000:.1f} pc ({frame['kpcPerPx'] * 1000 * 3.2616:.0f} ly) a pixel")
