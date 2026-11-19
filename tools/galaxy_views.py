"""Gives a galaxy of galaxy_fit.py what its photograph cannot: how it looks from above, finer
than the photograph's tilt allows, and how thick it is. From two pictures made for it (by an image
model, from the photograph: one face on, one edge on), and the photograph itself, which stays the
truth for what it shows.

A photograph of a tilted disc is squeezed along the disc's short axis: Andromeda, tilted by 70
degrees, has a pixel of it cover four times as much of the disc that way as the other way. What
lies within that, and how thick the disc is, a photograph does not say. So:

  face matrix  The face-on picture is turned, mirrored and scaled until, seen from the Earth, it
               lies as the photograph's galaxy does. It is kept as it is: the photograph only says
               where the disc lies, the picture what it looks like.
  side matrix  The edge-on picture, laid level and scaled to the disc's radius: at each radius,
               the light at each height, as drawn.
  lens         A few numbers that thicken the middle a little, as a convex lens.

Stars are then (viewer/src/fitted.rs): where on the disc from the face matrix, how high from the
side matrix's row at their radius, pushed apart by the lens.

Usage: galaxy_views.py FITTED.json NAME --face FILE --edge FILE [--size 1024]
Writes NAME_face.png next to FITTED.json (the disc from above: rows along its long axis, columns
along its short one; colours as light^(1/2.2)) and puts "faceMap", "sideMap" and "lens" into its
entry.
"""
import argparse
import json
import os

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates, rotate

args = argparse.ArgumentParser()
args.add_argument("fitted")
args.add_argument("name")
args.add_argument("--face", required=True)
args.add_argument("--edge", required=True)
args.add_argument("--size", type=int, default=1024)
args.add_argument("--radius-from-light", action="store_true", help="the disc's radius from where the photograph has light (scattered knots too), not from its bright body")
args.add_argument("--face-frame", type=float, help="the face-on picture is already laid on the disc (made from the photograph stretched face on): centred, rows along the long axis, columns along the short one, this many kpc from its middle to its edge; no search")
args.add_argument("--photo-angle", type=float, help="the face-on picture is the photograph stretched face on in its own orientation (the photograph's long axis where it is, only the short axis stretched back), centred, this many kpc from its middle to its side edge; no search")
args.add_argument("--min-likeness", type=float, default=0.0, help="refuse a face-on picture that, seen from the Earth, is less like the photograph than this (nothing written)")
args.add_argument("--reuse-turn", action="store_true", help="keep the turn and mirror found before, search only the scale")
args.add_argument("--data", default=os.path.expanduser("~/.gaiasky/data"))
args.add_argument("--preview")
a = args.parse_args()

doc = json.load(open(a.fitted))
base = os.path.dirname(os.path.abspath(a.fitted))
entry = next(g for g in doc["galaxies"] if a.name in g["names"] or g["name"] == a.name)
meta = entry["photoMetadata"]
photo_file = entry["photo"].replace("$data", a.data) if entry["photo"].startswith("$data") else os.path.join(base, entry["photo"])


def linear(file):
    return (np.asarray(Image.open(file).convert("RGB"), np.float32) / 255.0) ** 2.2


def without_sky(rgb):
    lum = rgb.mean(-1)
    edge = np.concatenate([lum[0], lum[-1], lum[:, 0], lum[:, -1]])
    return np.clip(rgb - np.median(edge) * 1.5, 0.0, None)


# ---- the photograph and where its pixels lie on the disc
photo = without_sky(linear(photo_file))
H, W = photo.shape[:2]
cx, cy = meta["middlePx"]
kpp = meta["kpcPerPx"]
(lx, ly), (sx, sy) = meta["longAxisPx"], meta["shortAxisPx"]
cos_i = meta["cosTilt"]
radius = meta["reachKpc"] / 1.15         # the disc's radius, as galaxy_fit.py measured it
if a.radius_from_light:
    # the disc as far out as the photograph has light on it, put back on the disc: galaxies whose
    # outer parts are scattered knots (M33) reach much further than their bright body
    _rgb = linear(photo_file)
    _l = _rgb.mean(-1)
    _yy, _xx = np.mgrid[0:_l.shape[0], 0:_l.shape[1]]
    _dx, _dy = _xx - meta["middlePx"][0], _yy - meta["middlePx"][1]
    _R = np.hypot((_dx * meta["longAxisPx"][0] + _dy * meta["longAxisPx"][1]) * meta["kpcPerPx"], (_dx * meta["shortAxisPx"][0] + _dy * meta["shortAxisPx"][1]) * meta["kpcPerPx"] / meta["cosTilt"])
    _lit = gaussian_filter(_l, 2) > 0.03 * np.quantile(_l, 0.999)
    measured = float(np.quantile(_R[_lit], 0.98))
    print(f"  the disc's radius from where the photograph has light: {measured:.2f} kpc (it was {radius:.2f})")
    radius = measured
# the face matrix reaches well beyond it, so that nothing of the picture is cut off: the galaxy's
# outline is where the picture has light, not a circle
reach = max(1.6 * radius, a.face_frame or 0.0, (a.photo_angle or 0.0) * 1.42)     # and the whole of a picture laid on the disc
N = a.size
mpx = 2 * reach / N                      # kpc a pixel of the map
g = (np.arange(N) + 0.5) * mpx - reach
X, Z = np.meshgrid(g, g, indexing="ij")  # rows along the long axis, columns along the short one
inside = np.ones_like(X, dtype=bool)
# where each place of the map is in the photograph
pcol = cx + (X * lx + Z * cos_i * sx) / kpp
prow = cy + (X * ly + Z * cos_i * sy) / kpp
photo_on_map = np.stack([map_coordinates(photo[..., c], [prow, pcol], order=1, cval=0.0) for c in range(3)], -1)

# ---- the face-on picture: its middle and radius
face = without_sky(linear(a.face))
fl = gaussian_filter(face.mean(-1), 4)
fcy, fcx = np.unravel_index(fl.argmax(), fl.shape)
ys, xs = np.nonzero(fl > 0.03 * fl.max())
f_radius = np.quantile(np.hypot(xs - fcx, ys - fcy), 0.995)


def face_on_map(turn_deg, mirror, scale, img=face):
    """The picture laid on the map: turned, mirrored, its radius scaled to the disc's."""
    if a.photo_angle:
        # in the photograph's own orientation: a place (x, z) of the disc is at x along the
        # photograph's long axis and z along its short one, kpc_f a pixel
        fh, fw = img.shape[:2]
        kf = 2.0 * a.photo_angle / fw
        col = (fw - 1) / 2.0 + (X * lx + Z * sx) / kf
        row = (fh - 1) / 2.0 + (X * ly + Z * sy) / kf
        return np.stack([map_coordinates(img[..., c], [row, col], order=1, cval=0.0) for c in range(3)], -1)
    t = np.radians(turn_deg)
    u, v = X, (-Z if mirror else Z)
    k = f_radius / radius * scale
    col = fcx + k * (u * np.cos(t) - v * np.sin(t))
    row = fcy + k * (u * np.sin(t) + v * np.cos(t))
    return np.stack([map_coordinates(img[..., c], [row, col], order=1, cval=0.0) for c in range(3)], -1)


def seen_from_earth(face_map, step=4):
    """The map as the photograph sees it, thin (on every step-th pixel of the photograph)."""
    rr, cc = np.mgrid[0:H:step, 0:W:step]
    dx, dy = cc - cx, rr - cy
    x = (dx * lx + dy * ly) * kpp
    z = (dx * sx + dy * sy) * kpp / cos_i
    return map_coordinates(face_map, [(x + reach) / mpx - 0.5, (z + reach) / mpx - 0.5], order=1, cval=0.0)


target = photo.mean(-1)[::4, ::4]
ok = target > 0.02 * target.max()


def likeness(face_map):
    s = seen_from_earth(face_map)
    m = ok | (s > 0.02 * s.max())
    return float(np.corrcoef(np.sqrt(s[m]), np.sqrt(target[m]))[0, 1])


if a.face_frame:
    # placed exactly: rows (the long axis) -> X, columns (the short axis) -> Z, as face_on_map
    # reads a picture turned by 90 degrees and mirrored, its middle the picture's middle
    fh, fw = face.shape[:2]
    fcy, fcx = (fh - 1) / 2.0, (fw - 1) / 2.0
    f_radius = radius * (fw / 2.0) / a.face_frame

# ---- turn, mirror, scale: the best likeness, coarse then fine (with --reuse-turn, the turn
# found before and only the scale)
fl_small = face.mean(-1)
best = (-2.0, 0.0, False, 1.0)
before = entry.get("faceMap", {})
reuse = (a.reuse_turn and "turnDeg" in before) or bool(a.face_frame) or bool(a.photo_angle)
if a.face_frame or a.photo_angle:
    before = {"turnDeg": 90.0, "mirror": True}
for mirror in ((bool(before.get("mirror")),) if reuse else (False, True)):
    for turn in ((float(before["turnDeg"]),) if reuse else np.arange(0, 360, 6)):
        # the picture's galaxy as large as the photograph's disc (the scale is measured, not
        # searched: a search drifts to sizes the galaxy does not have)
        for scale in (1.0,):
            m = face_on_map(turn, mirror, scale, img=face[..., :1].repeat(3, -1)).mean(-1) * inside
            like = likeness(m)
            if like > best[0]:
                best = (like, turn, mirror, scale)
_, turn, mirror, scale = best
for t in ((turn,) if reuse else np.arange(turn - 6, turn + 6.01, 1.0)):
    for sc in (scale,):
        m = face_on_map(t, mirror, sc, img=face[..., :1].repeat(3, -1)).mean(-1) * inside
        like = likeness(m)
        if like > best[0]:
            best = (like, t, mirror, sc)
fit_like, turn, mirror, scale = best
print(f"{entry['name']}: the face-on picture turned {turn:.1f} deg{', mirrored' if mirror else ''}, scaled {scale:.2f}: seen from the Earth like the photograph {fit_like:.3f}")
if fit_like < a.min_likeness:
    raise SystemExit(f"  refused: {fit_like:.3f} is under {a.min_likeness}; nothing written")

# ---- the face matrix: the picture itself, laid on the disc's axes (the photograph is not used
# for anything but where the disc lies and how it is turned)
picture = face_on_map(turn, mirror, scale) * inside[..., None]

# ---- the side matrix: the edge-on picture laid level, scaled to the disc's radius; for each of
# SIDE_R radii (both sides of the middle together) the light at SIDE_H heights, as drawn
SIDE_R, SIDE_H = 32, 64
edge = without_sky(linear(a.edge)).mean(-1)
ey, ex = np.nonzero(gaussian_filter(edge, 3) > 0.03 * gaussian_filter(edge, 3).max())
w_ = edge[ey, ex]
mx, my = np.average(ex, weights=w_), np.average(ey, weights=w_)
cov = np.cov(np.stack([ex - mx, ey - my]), aweights=w_)
angle = np.degrees(np.arctan2(*np.linalg.eigh(cov)[1][:, 1][::-1]))
level = rotate(edge, angle, reshape=False, order=1)
ls = gaussian_filter(level, 3)
lcy, lcx = np.unravel_index(ls.argmax(), ls.shape)
ys, xs = np.nonzero(ls > 0.03 * ls.max())
e_half = np.quantile(np.abs(xs - lcx), 0.995)
e_kpc = radius / e_half                       # kpc a pixel of the edge-on picture
height = 0.3 * radius                         # the matrix reaches this far above and below
rows_kpc = (np.arange(level.shape[0]) - lcy) * e_kpc
cols_kpc = np.abs(np.arange(level.shape[1]) - lcx) * e_kpc
side = np.zeros((SIDE_R, SIDE_H))
hb = np.linspace(-height, height, SIDE_H + 1)
for k in range(SIDE_R):
    cols = (cols_kpc >= k * radius / SIDE_R) & (cols_kpc < (k + 1) * radius / SIDE_R)
    prof = level[:, cols].sum(1)
    side[k] = np.histogram(rows_kpc, bins=hb, weights=prof)[0]
side = side / np.maximum(side.max(1, keepdims=True), 1e-12)
side_bytes = np.round(side * 255).astype(int)
widths = [float(np.sqrt(np.sum(r * (0.5 * (hb[1:] + hb[:-1])) ** 2) / max(r.sum(), 1e-12))) for r in side]
print("  side: rms height (kpc) " + ", ".join(f"{(k + 0.5) * radius / SIDE_R:.1f}: {w:.2f}" for k, w in list(enumerate(widths))[::4]))

# ---- written
name = entry["name"].replace(" ", "_").replace("'", "")
top = np.quantile(picture.mean(-1)[inside], 0.9995)
Image.fromarray((np.clip(picture / top, 0, 1) ** (1 / 2.2) * 255 + 0.5).astype(np.uint8)).save(os.path.join(base, f"{name}_face.png"))
entry["faceMap"] = {"file": f"{name}_face.png", "halfKpc": float(reach), "size": N,
                    "about": "the face matrix: the disc from above, rows along its long axis, columns along its short one; colours as light^(1/2.2). The face-on picture as it is, turned onto the disc's axes by tools/galaxy_views.py",
                    "facePicture": os.path.relpath(os.path.abspath(a.face), base), "turnDeg": float(turn), "mirror": bool(mirror), "scale": float(scale),
                    "fromPhotograph": bool(a.photo_angle or a.face_frame), "likePhotograph": float(fit_like)}
entry["sideMap"] = {"radii": SIDE_R, "heights": SIDE_H, "radiusKpc": float(radius), "heightKpc": float(height),
                    "edgePicture": os.path.relpath(os.path.abspath(a.edge), base),
                    "about": "the side matrix: for each of 'radii' radii from the middle out to radiusKpc, the light at 'heights' heights from -heightKpc to +heightKpc (0-255, each row's brightest 255), from the edge-on picture as drawn",
                    "values": side_bytes.ravel().tolist()}
# the lens: the disc thicker towards the middle as a convex lens, and some light of its own there
# (a picture's middle is burnt out), in proportion to the disc
entry["lens"] = {"radiusKpc": float(0.3 * radius), "strength": 0.4, "extraLight": 0.0, "extraRadiusKpc": float(0.15 * radius), "extraFlat": 0.2,
                 "about": "within radiusKpc the heights are pushed apart by up to 1 + strength in the middle (y + strength dome h tanh(y / h), dome = sqrt(1 - (R / radiusKpc)^2), h = (0.12 dome + 0.02) of the radius); extraLight of the stars in a flattened Gaussian of extraRadiusKpc"}
entry.pop("thickness", None)
json.dump(doc, open(a.fitted, "w"), indent=2)
print("wrote", os.path.join(base, f"{name}_face.png"), "and the entry in", a.fitted)
if a.preview:
    s = lambda v: (np.clip(v / np.quantile(v, 0.999), 0, 1) ** (1 / 2.2) * 255).astype(np.uint8)
    Image.fromarray(np.concatenate([s(picture), s(photo_on_map * inside[..., None])], 1)).save(f"{a.preview}/{name}_face_steps.png")
