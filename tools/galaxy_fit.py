"""Fits a spiral galaxy's formula to its photograph: a few dozen numbers instead of voxel grids.

The star formula of the viewer (viewer/src/procedural.rs, ParamGalaxy) makes the stars of a galaxy
from these numbers alone, so that a catalog of many galaxies costs a few hundred bytes each. This
tool finds the numbers from a photograph by itself, pixel by pixel, with no one to look at it:

  1. the ellipse of the disc (as photo_galaxy.py): its middle, the direction of its long axis and
     how flat it looks, which gives how it is tilted; with the place of the galaxy in the sky
     (from the catalog) that gives where its axes point in space;
  2. the disc as from above: the ellipse stretched back into a circle;
  3. the light along the radius, as a round bulge (a Sersic profile: its share, radius and
     steepness) and an exponential disc (its scale and where it ends);
  4. the spiral arms, from what is left once the light is divided by that: a Fourier transform in
     the angle and the logarithm of the radius finds how many arms there are (1 to 4), how
     tightly they wind (the pitch), where they start and how strong they are; their sharpness and
     the dust lanes along their inner edges by trying a few;
  5. the colours of the middle, of the disc and of the arms.

It then draws the formula as seen from the Earth and says how much it looks like the photograph
(the correlation of the square roots of the light, pixel by pixel; 1 is the same).

Usage: galaxy_fit.py OUT.json NAME [NAME ...] [--data DIR] [--radius-kpc NAME=R ...] [--preview DIR]
A NAME is one of the names of a galaxy in the data's catalog-nbg (M33, M81, Andromeda Galaxy, ...),
which gives its picture, where it is and how big. The picture's width is its size in the catalog,
or, with --radius-kpc, its disc reaches R kpc (99.5% of it, as photo_galaxy.py measures).
"""
import argparse
import glob
import json
import os

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares

args = argparse.ArgumentParser()
args.add_argument("out")
args.add_argument("names", nargs="+")
args.add_argument("--data", default=os.path.expanduser("~/.gaiasky/data"))
args.add_argument("--radius-kpc", nargs="*", default=[])
args.add_argument("--photo", nargs="*", default=[], help="NAME=FILE: another photograph than the catalog's (an upscaled one), the same view")
args.add_argument("--abs-mag", nargs="*", default=[], help="NAME=M where the catalog has no magnitude")
args.add_argument("--core-only", nargs="*", default=[], help="NAME: measure the disc on its body around the middle only (a companion galaxy is joined to it)")
args.add_argument("--radec", nargs="*", default=[], help="NAME=RA,DEC (degrees): where the catalog has it wrong (catalog-nbg: M65's right ascension in hours, M58's and M61's off by degrees)")
args.add_argument("--tilt", nargs="*", default=[], help="NAME=DEG: the tilt where the light's outline does not show it (a disc inside a big round halo: the Sombrero)")
args.add_argument("--preview")
a = args.parse_args()
radius_given = dict((k, float(v)) for k, v in (r.split("=") for r in a.radius_kpc))
photo_given = dict(r.split("=", 1) for r in a.photo)
core_only = set(a.core_only)
radec_given = dict((k, tuple(float(x) for x in v.split(","))) for k, v in (r.split("=") for r in a.radec))
tilt_given = dict((k, float(v)) for k, v in (r.split("=") for r in a.tilt))
mag_given = dict((k, float(v)) for k, v in (r.split("=") for r in a.abs_mag))

N = 256                       # pixels across the galaxy as from above
HEIGHT_PER_RADIUS = 0.005     # the disc's thickness over its radius (not in a photograph)
BULGE_FLAT = 0.7              # how round the bulge is, short over long (not in a photograph)


def smoothstep(lo, hi, v):
    t = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def rot(which, deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return {"y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[which]


TO_EQ = rot("y", 90 + 192.85948) @ rot("z", 90 - 27.12825) @ rot("y", -32.93192)    # galactic -> equatorial


def catalog(data):
    """Every galaxy of catalog-nbg's JSON files with a place: names -> the object."""
    found = {}
    def walk(o):
        if isinstance(o, dict):
            eq = (o.get("coordinates") or {}).get("equatorial")
            if o.get("names") and eq:
                for n in o["names"]:
                    found[n.replace(" ", "").lower()] = o
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    for f in glob.glob(f"{data}/catalog-nbg/*.json"):
        walk(json.load(open(f)))
    return found


def picture_of(o, data):
    tex = ((o.get("model") or {}).get("material") or {}).get("diffuse")
    if tex:
        return tex.replace("$data", data)
    for n in o["names"]:
        f = f"{data}/catalog-nbg/tex/extragal/{n.replace(' ', '')}.jpg"
        if os.path.exists(f):
            return f
    raise SystemExit(f"no picture for {o['names'][0]}")


def sersic_b(n):
    return 2.0 * n - 1.0 / 3.0 + 0.009876 / n


def fit(name, o):
    photo_file = next((photo_given[n] for n in [name] + o["names"] if n in photo_given), None) or picture_of(o, a.data)
    ra, dec, dist_pc = o["coordinates"]["equatorial"]
    radec = next((radec_given[n] for n in [name] + o["names"] if n in radec_given), None)
    if radec is not None:
        ra, dec = radec
    photo = (np.asarray(Image.open(photo_file).convert("RGB"), np.float32) / 255.0) ** 2.2
    H, W = photo.shape[:2]
    lum = photo.mean(axis=-1)

    # ---- 1. the ellipse
    smooth = gaussian_filter(lum, 6 * W / 1024)
    cy, cx = np.unravel_index(smooth.argmax(), smooth.shape)
    outline = smooth > 0.02 * smooth.max()
    if name in core_only or any(n in core_only for n in o["names"]):
        # only the body around the brightest place: a companion joined to it by a faint bridge
        # (M51's NGC 5195) is not part of the disc. The body: brighter than 12% of the middle,
        # in one piece with it
        from scipy.ndimage import label
        pieces, _ = label(smooth > 0.12 * smooth.max())
        body = pieces == pieces[cy, cx]
        grown = gaussian_filter(body.astype(np.float32), 25 * W / 1024) > 0.05
        outline = outline & grown
    ys, xs = np.nonzero(outline)
    values, vectors = np.linalg.eigh(np.cov(np.stack([xs - cx, ys - cy])))
    long_axis = vectors[:, 1] * np.sign(vectors[0, 1] or 1.0)
    short_axis = np.array([-long_axis[1], long_axis[0]])
    along = (xs - cx) * long_axis[0] + (ys - cy) * long_axis[1]
    reach_px, flat = np.quantile(np.abs(along), 0.995), float(np.sqrt(values[0] / values[1]))
    if name in radius_given or any(n in radius_given for n in o["names"]):
        radius_kpc = radius_given.get(name) or next(radius_given[n] for n in o["names"] if n in radius_given)
        kpc_per_px = radius_kpc / reach_px
    else:
        kpc_per_px = o["sizePc"] / 1000.0 / W
        radius_kpc = reach_px * kpc_per_px
    # a disc of some thickness never looks flatter than it: the tilt from the axis ratio
    q0 = 0.12
    cos_i = np.sqrt(max(flat * flat - q0 * q0, 0.0) / (1 - q0 * q0))
    tilt = next((tilt_given[n] for n in [name] + o["names"] if n in tilt_given), None)
    if tilt is not None:
        cos_i = float(np.cos(np.radians(tilt)))
    # where its axes point (as photo_galaxy.py; the near side is not in a photograph: the north-west)
    r_, d_ = np.radians(ra), np.radians(dec)
    sight = np.array([np.cos(d_) * np.sin(r_), np.sin(d_), np.cos(d_) * np.cos(r_)])
    north = np.array([0.0, 1.0, 0.0]) - sight[1] * sight
    north /= np.linalg.norm(north)
    east = np.cross(north, sight)
    sky = lambda v: -v[0] * east - v[1] * north
    long_sky, short_sky = sky(long_axis), sky(short_axis)
    near = north - east
    leans_away = 1.0 if short_sky @ near < 0 else -1.0
    short_in_space = cos_i * short_sky + leans_away * np.sqrt(1 - cos_i * cos_i) * sight
    pole = np.cross(long_sky, short_in_space)
    pole /= np.linalg.norm(pole)
    axes = [TO_EQ.T @ w for w in (long_sky, pole, short_in_space)]
    centre_gal = TO_EQ.T @ (sight * dist_pc)

    # ---- 2. as from above: x along the long axis, z along the short one (kpc)
    half = radius_kpc * 1.15
    g = (np.arange(N) + 0.5) / N * 2 * half - half
    X, Z = np.meshgrid(g, g, indexing="ij")
    col = cx + (X * long_axis[0] + Z * cos_i * short_axis[0]) / kpc_per_px
    row = cy + (X * long_axis[1] + Z * cos_i * short_axis[1]) / kpc_per_px
    face = np.stack([map_coordinates(gaussian_filter(photo[..., c], 0.7), [row, col], order=1, cval=0.0) for c in range(3)], -1) * cos_i
    face_l = face.mean(-1)
    R = np.hypot(X, Z)
    TH = np.arctan2(Z, X)

    # ---- 3. the light along the radius: bulge (Sersic) + disc (exponential)
    edges = np.linspace(0, radius_kpc, 60)
    rr = 0.5 * (edges[1:] + edges[:-1])
    prof = np.array([np.median(face_l[(R >= edges[k]) & (R < edges[k + 1])]) if ((R >= edges[k]) & (R < edges[k + 1])).any() else 0.0 for k in range(len(rr))])
    keep = prof > prof.max() * 1e-3
    def model_prof(p, r):
        ib, re, n, idisc, h = np.exp(p[0]), np.exp(p[1]), p[2], np.exp(p[3]), np.exp(p[4])
        return ib * np.exp(-sersic_b(n) * ((r / re) ** (1 / n) - 1)) + idisc * np.exp(-r / h)
    p0 = [np.log(prof.max() * 0.3), np.log(radius_kpc * 0.08), 2.0, np.log(prof.max() * 0.5), np.log(radius_kpc * 0.25)]
    # on the square root of the light: the bright middle weighs as much as the faint edge
    res = least_squares(lambda p: np.sqrt(model_prof(p, rr[keep])) - np.sqrt(prof[keep]), p0,
                        bounds=([-20, np.log(radius_kpc * 0.005), 0.5, -20, np.log(radius_kpc * 0.03)], [5, np.log(radius_kpc * 0.5), 6.0, 5, np.log(radius_kpc * 2)]))
    ib, re, n, idisc, h = np.exp(res.x[0]), np.exp(res.x[1]), res.x[2], np.exp(res.x[3]), np.exp(res.x[4])
    # where the disc ends: its light falls below 3% of where it is at one scale length
    disc_light = idisc * np.exp(-rr / h)
    edge = float(min(radius_kpc * 1.1, rr[np.nonzero(prof > 0.03 * idisc * np.exp(-1))[0].max()] if (prof > 0.03 * idisc * np.exp(-1)).any() else radius_kpc))
    # the shares of the light, as seen from above (both integrated over the plane)
    rfine = np.linspace(0, edge, 2000)
    dr = rfine[1] - rfine[0]
    bulge_l = np.sum(ib * np.exp(-sersic_b(n) * ((rfine / re) ** (1 / n) - 1)) * 2 * np.pi * rfine) * dr
    disc_l = np.sum(idisc * np.exp(-rfine / h) * smoothstep(edge, edge * 0.85, rfine) * 2 * np.pi * rfine) * dr
    bulge_share = float(bulge_l / (bulge_l + disc_l))
    axisym = model_prof(res.x, R) * smoothstep(edge, edge * 0.85, R)

    # ---- 4. the arms: the light over the round model, in angle and log radius
    ring = (R > max(2.0 * re, 0.12 * edge)) & (R < 0.9 * edge)
    ratio = np.where(ring, np.log((face_l + 1e-4) / (axisym + 1e-4)), 0.0)
    ratio = np.where(ring, ratio - ratio[ring].mean(), 0.0)
    w_ = ring * R ** 0   # each pixel the same
    u = np.log(np.maximum(R, 1e-3))
    best = (0.0, 2, 20.0, 0.0)
    for m in (2, 3, 4):
        for pitch in np.arange(5.0, 45.0, 0.5):
            k = m / np.tan(np.radians(pitch))
            for sense in (1.0, -1.0):
                c = np.sum(ratio * w_ * np.exp(-1j * (m * TH - sense * k * u))) / max(np.sum(w_), 1)
                # a pattern of more arms is less likely by chance: the weight of one arm
                score = abs(c)
                if score > best[0]:
                    best = (score, m, pitch * sense, float(np.angle(c)))
    amp, m, pitch, phase = min(2.0 * best[0], 1.0), best[1], best[2], best[3]
    k = m / np.tan(np.radians(abs(pitch))) * np.sign(pitch)
    psi = m * TH - k * u - phase
    start = float(max(2.0 * re, 0.12 * edge))

    def arms_of(amp_, sharp, dust, offset, psi_=psi, r_=R):
        s = ((1 + np.cos(psi_)) / 2) ** sharp
        mean_s = np.exp(np.sum(np.log(np.arange(sharp + 1, 2 * sharp + 1) / np.arange(1, sharp + 1)))) / 4 ** sharp if sharp == int(sharp) else s.mean()
        on = smoothstep(start * 0.6, start * 1.4, r_)
        arms = 1 + amp_ * on * (s / mean_s - 1)
        lane = 1 - dust * on * ((1 + np.cos(psi_ - offset)) / 2) ** 8
        return np.maximum(arms, 0.0) * lane
    target = np.where(R < edge, face_l, 0.0)
    best2 = None
    for sharp in (1, 2, 4, 8):
        for amp_s in (0.25, 0.5, 0.75, 1.0) if sharp > 1 else (0.5, 1.0, 1.5, 2.0):
            for dust in (0.0, 0.3, 0.6):
                for offset in ((0.0,) if dust == 0 else (-0.8, -0.4, 0.4, 0.8)):
                    model = axisym * arms_of(min(amp * amp_s * 2, 3.0), sharp, dust, offset)
                    model *= target.sum() / max(model.sum(), 1e-12)
                    e = np.mean((np.sqrt(model) - np.sqrt(target)) ** 2)
                    if best2 is None or e < best2[0]:
                        best2 = (e, min(amp * amp_s * 2, 3.0), sharp, dust, offset)
    _, strength, sharp, dust, dust_offset = best2

    # ---- 5. what is left, coarsely: a 32 x 32 map of the light over the formula, as from above
    # (lopsided discs, big dark lanes, clouds of young stars), in steps of 1/16 of e
    model_face = axisym * arms_of(strength, sharp, dust, dust_offset)
    model_face *= target.sum() / max(model_face.sum(), 1e-12)
    D = 32
    inside = R < edge
    num = (np.where(inside, np.log((face_l + 1e-4) / (model_face + 1e-4)), 0.0)).reshape(D, N // D, D, N // D).sum((1, 3))
    cnt = inside.reshape(D, N // D, D, N // D).sum((1, 3))
    detail = np.clip(np.where(cnt > 0, num / np.maximum(cnt, 1), 0.0), -2.5, 2.5)
    detail = np.round(detail * 16).astype(int)

    def detail_at(x_kpc, z_kpc):
        f = map_coordinates(detail / 16.0, [(x_kpc + half) / (2 * half) * D - 0.5, (z_kpc + half) / (2 * half) * D - 0.5], order=1, mode="nearest")
        return np.exp(f)

    # ---- 6. colours (linear)
    def mean_colour(mask):
        return [float(v) for v in (face[mask].mean(0) / max(face[mask].mean(), 1e-9))] if mask.any() else [1.0, 1.0, 1.0]
    on_arm = arms_of(strength, sharp, 0.0, 0.0) > 1.3
    colours = {"middle": mean_colour(R < re), "disc": mean_colour(ring & ~on_arm), "arms": mean_colour(ring & on_arm)}

    absmag = next((mag_given[n] for n in [name] + o["names"] if n in mag_given), None)
    if absmag is None:
        absmag = o.get("absMag")
    if absmag is None and o.get("appMag"):
        absmag = o["appMag"] - 5 * np.log10(dist_pc / 10)
    absmag = float(absmag if absmag is not None else -20.0)
    # as many systems as its light says, the Milky Way's 2.7e11 for absolute magnitude -20.8
    systems = 2.7e11 * 10 ** (-0.4 * (absmag + 20.8))

    entry = {
        "name": o["names"][0],
        "names": o["names"],
        # the catalog's photograph as $data/...; another one relative to the output file
        "photo": photo_file.replace(a.data, "$data") if photo_file.startswith(a.data) else os.path.relpath(os.path.abspath(photo_file), os.path.dirname(os.path.abspath(a.out))),
        "photoMetadata": {"sizePx": [W, H], "middlePx": [int(cx), int(cy)], "kpcPerPx": float(kpc_per_px),
                          "longAxisDeg": float(np.degrees(np.arctan2(-long_axis[0], -long_axis[1])) % 180), "axisRatio": flat,
                          "longAxisPx": [float(long_axis[0]), float(long_axis[1])], "shortAxisPx": [float(short_axis[0]), float(short_axis[1])],
                          "cosTilt": float(cos_i), "reachKpc": float(half),
                          "about": "pixel (col, row) of a place (x, z) of the disc, kpc along its long and short axes: middlePx + (x longAxisPx + z cosTilt shortAxisPx) / kpcPerPx; the photograph's light there is how many stars there are"},
        "centrePc": [float(v) for v in centre_gal],
        "axes": [[float(v) for v in w] for w in axes],
        "tiltDeg": float(np.degrees(np.arccos(cos_i))),
        "absMag": absmag,
        "systems": float(systems),
        "bulge": {"share": bulge_share, "radiusKpc": float(re), "sersic": float(n), "flattening": BULGE_FLAT},
        "disc": {"scaleKpc": float(h), "heightKpc": float(HEIGHT_PER_RADIUS * radius_kpc), "edgeKpc": edge},
        "arms": {"count": int(m), "pitchDeg": float(pitch), "phaseDeg": float(np.degrees(phase)), "strength": float(strength), "sharpness": int(sharp), "startKpc": start},
        "dust": {"strength": float(dust), "offsetDeg": float(np.degrees(dust_offset))},
        "colours": colours,
        "detail": {"size": D, "halfKpc": float(half), "about": "rows along the long axis, columns along the short one; the natural log of the light over the formula, in sixteenths", "values": detail.ravel().tolist()},
    }

    # ---- how like the photograph: the formula as from above, seen again from the Earth
    model_face = axisym * arms_of(strength, sharp, dust, dust_offset) * detail_at(X, Z)
    rr_px, cc_px = np.mgrid[0:H, 0:W]
    dx, dy = cc_px - cx, rr_px - cy
    xk = (dx * long_axis[0] + dy * long_axis[1]) * kpc_per_px
    zk = (dx * short_axis[0] + dy * short_axis[1]) * kpc_per_px / max(cos_i, 0.05)
    seen = map_coordinates(model_face, [(xk + half) / (2 * half) * N - 0.5, (zk + half) / (2 * half) * N - 0.5], order=1, cval=0.0)
    mask = (seen > 0.01 * seen.max()) | (lum > 0.02 * lum.max())
    like = float(np.corrcoef(np.sqrt(seen[mask]), np.sqrt(lum[mask]))[0, 1])
    entry["fit"] = {"likePhotograph": like}
    print(f"{entry['name']}: tilted {entry['tiltDeg']:.0f} deg, bulge {bulge_share:.2f} of the light (r {re:.2f} kpc, n {n:.1f}), disc {h:.2f} kpc to {edge:.1f} kpc, "
          f"{m} arms pitch {pitch:.1f} deg strength {strength:.2f} sharp {sharp}, dust {dust:.1f}; like the photograph {like:.3f}")
    if a.preview:
        os.makedirs(a.preview, exist_ok=True)
        s = lambda v: (np.clip(v / np.quantile(v, 0.999), 0, 1) ** 0.5 * 255).astype(np.uint8)
        Image.fromarray(np.concatenate([s(lum), s(seen)], 1)).save(f"{a.preview}/{entry['name'].replace(' ', '_')}_fit.png")
    return entry


cat = catalog(a.data)
out = []
for name in a.names:
    o = cat.get(name.replace(" ", "").lower())
    if o is None:
        raise SystemExit(f"no galaxy called {name} in {a.data}/catalog-nbg")
    out.append(fit(name, o))
json.dump({"about": "Galaxies as formulas, fitted to their photographs by tools/galaxy_fit.py (see there for what each number is). Lengths in kpc along the galaxy's own axes (axes: the long axis, the pole, the short axis, galactic); centrePc galactic, parsecs. photoMetadata: where the photograph has the galaxy.", "galaxies": out}, open(a.out, "w"), indent=2)
print("wrote", a.out)
