"""Galaxies without a photograph, from two pictures imagined for them and the catalog's numbers.

For the galaxies of Gaia Sky's catalog-nbg that have no photograph, an image model made a face-on
and an edge-on picture from what the catalog says (name, size, brightness): the famous ones as
they are, the others as galaxies of that size and brightness. This tool turns them into entries of
the same kind as galaxy_views.py makes (viewer/src/fitted.rs reads both):

  where       the catalog's right ascension, declination and distance;
  how turned  there is no photograph to say: a direction of its pole drawn from its name, so that
              it is always the same;
  how big     the catalog's diameter;
  how many    stars as its absolute magnitude says (the Milky Way's 2.7e11 systems at -20.8);
  face matrix the face-on picture, its galaxy's radius laid on the catalog's;
  side matrix the edge-on picture laid level, as galaxy_views.py does;
  lens        as galaxy_views.py.

Usage: galaxy_imagine.py OUT.json GALAXIES.json PICTURES_DIR [--size 512]
GALAXIES.json: a list of {name, ra, dec, distancePc, absMag, diameterKpc}; PICTURES_DIR has
NAME_face.png and NAME_edge.png for each. Writes OUT.json and the face matrices next to it in a
folder named as OUT without its extension.
"""
import argparse
import hashlib
import json
import os

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates, rotate

args = argparse.ArgumentParser()
args.add_argument("out")
args.add_argument("galaxies")
args.add_argument("pictures")
args.add_argument("--size", type=int, default=512)
a = args.parse_args()

SIDE_R, SIDE_H = 32, 64


def rot(which, deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return {"y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[which]


TO_EQ = rot("y", 90 + 192.85948) @ rot("z", 90 - 27.12825) @ rot("y", -32.93192)    # galactic -> equatorial


def linear(file):
    rgb = (np.asarray(Image.open(file).convert("RGB"), np.float32) / 255.0) ** 2.2
    lum = rgb.mean(-1)
    sky = np.median(np.concatenate([lum[0], lum[-1], lum[:, 0], lum[:, -1]])) * 1.5
    return np.clip(rgb - sky, 0, None)


def seeded(name):
    return np.random.default_rng(int(hashlib.sha256(name.encode()).hexdigest()[:12], 16))


def entry_for(g, folder, out_dir):
    name = g["name"]
    key = name.replace(" ", "_")
    radius = g["diameterKpc"] / 2.0
    # ---- where and how turned
    ra, dec = np.radians(g["ra"]), np.radians(g["dec"])
    sight = np.array([np.cos(dec) * np.sin(ra), np.sin(dec), np.cos(dec) * np.cos(ra)])
    centre_gal = TO_EQ.T @ (sight * g["distancePc"])
    rng = seeded(name)
    pole = rng.normal(size=3)
    pole /= np.linalg.norm(pole)
    long_ = np.cross(pole, [0.0, 0.0, 1.0] if abs(pole[2]) < 0.9 else [1.0, 0.0, 0.0])
    long_ /= np.linalg.norm(long_)
    short = np.cross(long_, pole)
    # ---- the face matrix: the picture's galaxy as wide as the catalog's
    face = linear(f"{folder}/{key}_face.png")
    fl = gaussian_filter(face.mean(-1), 4)
    fcy, fcx = np.unravel_index(fl.argmax(), fl.shape)
    ys, xs = np.nonzero(fl > 0.03 * fl.max())
    f_radius = np.quantile(np.hypot(xs - fcx, ys - fcy), 0.995)
    half = 1.25 * radius
    n = a.size
    gr = (np.arange(n) + 0.5) / n * 2 * half - half
    X, Z = np.meshgrid(gr, gr, indexing="ij")
    k = f_radius / radius
    m = np.stack([map_coordinates(face[..., c], [fcy + k * X, fcx + k * Z], order=1, cval=0.0) for c in range(3)], -1)
    top = np.quantile(m.mean(-1), 0.9995)
    face_file = f"{key}_face.png"
    Image.fromarray((np.clip(m / max(top, 1e-9), 0, 1) ** (1 / 2.2) * 255 + 0.5).astype(np.uint8)).save(os.path.join(out_dir, face_file))
    # ---- the side matrix (as galaxy_views.py)
    edge = linear(f"{folder}/{key}_edge.png").mean(-1)
    es = gaussian_filter(edge, 3)
    ey, ex = np.nonzero(es > 0.03 * es.max())
    w_ = edge[ey, ex]
    mx, my = np.average(ex, weights=w_), np.average(ey, weights=w_)
    vec = np.linalg.eigh(np.cov(np.stack([ex - mx, ey - my]), aweights=w_))[1][:, 1]
    level = rotate(edge, np.degrees(np.arctan2(vec[1], vec[0])), reshape=False, order=1)
    ls = gaussian_filter(level, 3)
    lcy, lcx = np.unravel_index(ls.argmax(), ls.shape)
    ys, xs = np.nonzero(ls > 0.03 * ls.max())
    e_kpc = radius / max(np.quantile(np.abs(xs - lcx), 0.995), 1.0)
    height = 0.3 * radius
    rows_kpc = (np.arange(level.shape[0]) - lcy) * e_kpc
    cols_kpc = np.abs(np.arange(level.shape[1]) - lcx) * e_kpc
    hb = np.linspace(-height, height, SIDE_H + 1)
    side = np.zeros((SIDE_R, SIDE_H))
    for r in range(SIDE_R):
        cols = (cols_kpc >= r * radius / SIDE_R) & (cols_kpc < (r + 1) * radius / SIDE_R)
        side[r] = np.histogram(rows_kpc, bins=hb, weights=level[:, cols].sum(1))[0]
    side = side / np.maximum(side.max(1, keepdims=True), 1e-12)
    systems = 2.7e11 * 10 ** (-0.4 * (g["absMag"] + 20.8))
    return {
        "name": name, "names": [name], "photo": "", "imagined": True,
        "about": "imagined: no photograph; face-on and edge-on pictures made by an image model from the catalog's numbers (tools/galaxy_imagine.py)",
        "catalog": {k2: g[k2] for k2 in ("ra", "dec", "distancePc", "absMag", "diameterKpc")},
        "centrePc": [float(v) for v in centre_gal],
        "axes": [[float(v) for v in TO_EQ.T @ w] for w in (long_, pole, short)],
        "absMag": float(g["absMag"]), "systems": float(systems),
        # the formula's numbers (for a galaxy without its matrices): plain
        "bulge": {"share": 0.0, "radiusKpc": 0.1 * radius, "sersic": 1.0, "flattening": 0.7},
        "disc": {"scaleKpc": radius / 3.0, "heightKpc": 0.015 * radius, "edgeKpc": radius},
        "arms": {"count": 2, "pitchDeg": 20.0, "phaseDeg": 0.0, "strength": 0.0, "sharpness": 2, "startKpc": 0.2 * radius},
        "dust": {"strength": 0.0, "offsetDeg": 0.0},
        "colours": {"middle": [1, 1, 1], "disc": [1, 1, 1], "arms": [1, 1, 1]},
        "faceMap": {"file": os.path.join(os.path.basename(out_dir), face_file), "halfKpc": float(half), "size": n,
                    "about": "the face-on picture, its galaxy's radius laid on the catalog's", "facePicture": f"{key}_face.png"},
        "sideMap": {"radii": SIDE_R, "heights": SIDE_H, "radiusKpc": float(radius), "heightKpc": float(height), "edgePicture": f"{key}_edge.png",
                    "values": np.round(side * 255).astype(int).ravel().tolist()},
        "lens": {"radiusKpc": float(0.3 * radius), "strength": 0.4, "extraLight": 0.12, "extraRadiusKpc": float(0.15 * radius), "extraFlat": 0.2},
    }


galaxies = json.load(open(a.galaxies))
out_dir = os.path.splitext(a.out)[0]
os.makedirs(out_dir, exist_ok=True)
entries = []
for g in galaxies:
    key = g["name"].replace(" ", "_")
    if not (os.path.exists(f"{a.pictures}/{key}_face.png") and os.path.exists(f"{a.pictures}/{key}_edge.png")):
        print("no pictures for", g["name"])
        continue
    entries.append(entry_for(g, a.pictures, out_dir))
json.dump({"about": "Galaxies without a photograph, imagined from the catalog (tools/galaxy_imagine.py); the same kind of entries as fitted.json.", "galaxies": entries}, open(a.out, "w"))
print(f"wrote {len(entries)} galaxies to {a.out}")
