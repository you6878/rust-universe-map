"""A page to look at what stars a galaxy's three views make at a place (tools/galaxy_probe.html):
click one of the three pictures, or type the place in the game's coordinates (the galaxy's own:
X its long axis, Y its pole, Z its short axis, parsecs from its middle), and the stars of a cube
around it are drawn. The same density as viewer/src/photo_volume.rs, on the same frame
(tools/galaxy_xyz.py --photo: seen from the Earth like the photograph).

Usage: galaxy_probe.py FRAME.json OUT.html [PAGE.html]   (FRAME.json and its three pictures, from galaxy_xyz.py;
PAGE.html the page to fill, galaxy_probe.html by default, galaxy_cube.html for the voxel cube)
"""
import base64
import json
import os
import sys

frame_file, out = sys.argv[1], sys.argv[2]
meta = json.load(open(frame_file))
base = os.path.dirname(frame_file)


def data(view):
    with open(os.path.join(base, meta["views"][view]), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()


template = sys.argv[3] if len(sys.argv) > 3 else "galaxy_probe.html"
def light16(view):
    """the 16-bit light map (log10 of the light) as an RGB PNG: high byte red, low byte green"""
    import io
    import numpy as np
    from PIL import Image
    files = meta.get("light", {}).get("files", {})
    if view not in files:
        return ""
    v = np.asarray(Image.open(os.path.join(base, files[view]))).astype(np.uint32)
    rgb = np.stack([(v >> 8).astype(np.uint8), (v & 255).astype(np.uint8), np.zeros_like(v, np.uint8)], -1)
    buf = io.BytesIO()
    Image.fromarray(rgb).save(buf, "PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


FLANK, PEAK_MOST = 6, 2.0


def fill(p):
    """an edge view's column, its dust lane filled as viewer/src/photo_volume.rs fill"""
    import numpy as np
    n = len(p)
    flat = np.maximum(p, np.minimum(np.maximum.accumulate(p), np.maximum.accumulate(p[::-1])[::-1]))
    k = 0
    while k < n:
        if flat[k] <= p[k]:
            k += 1
            continue
        a = k
        while k < n and flat[k] > p[k]:
            k += 1
        b = k - 1
        level = flat[a]
        mid = (a + b) / 2
        sides = [j for j in list(range(max(0, a - FLANK), a)) + list(range(b + 1, min(n, b + 1 + FLANK))) if p[j] > 0]
        c = None
        if len(sides) >= 5:
            z = np.array(sides, float) - mid
            A = np.stack([np.ones_like(z), z, z * z], 1)
            sol, *_ = np.linalg.lstsq(A, np.log(p[sides].astype(float)), rcond=None)
            top = -sol[1] / (2 * sol[2]) if sol[2] != 0 else 0
            if sol[2] < 0 and z[0] < top < z[-1]:
                c = sol
        for j in range(a, b + 1):
            z = j - mid
            p[j] = min(level * PEAK_MOST, max(level, np.exp(c[0] + c[1] * z + c[2] * z * z))) if c is not None else level
    return p


def coarse(factor=8):
    """the density averaged over blocks of factor^3 voxels (stars a pc^3), 16-bit log10 (0 none),
    and the top view's colour over the same blocks; for the far voxels of galaxy_fly.html"""
    import io
    import numpy as np
    from PIL import Image
    lo, hi = meta["light"]["log10Range"]
    files = meta["light"]["files"]
    def light(view):
        v = np.asarray(Image.open(os.path.join(base, files[view]))).astype(np.float64)
        return np.where(v > 0, 10 ** (lo + (v - 1) / 65534 * (hi - lo)), 0.0)
    top, side_u, side_w = light("y"), light("x"), light("z")
    n = top.shape[0]
    FU = np.stack([fill(side_u[:, c].copy()) for c in range(n)])      # [col, h]
    FW = np.stack([fill(side_w[:, r].copy()) for r in range(n)])      # [row, h]
    norm = FW @ FU.T                                                  # [row, col]
    pc = meta["frame"]["kpcPerPx"] * 1000
    A = np.where(norm > 0, meta["stars"] * top / top.sum() / np.where(norm > 0, norm, 1), 0) / pc ** 3
    m = n // factor
    out = np.zeros((m, m, m))                                         # [h block, row block, col block]
    for h in range(n):
        d = FW[:, h, None] * A * FU[None, :, h]
        out[h // factor] += d.reshape(m, factor, m, factor).sum((1, 3))
    out /= factor ** 3
    print(f"coarse {m}^3: stars {out.sum() * (pc * factor) ** 3:.3e} (of {meta['stars']:.3e}), most {out.max():.3g} a pc^3")
    pos = out[out > 0]
    clo, chi = np.log10(pos.min()), np.log10(pos.max())
    enc = np.where(out > 0, np.round((np.log10(np.maximum(out, 1e-300)) - clo) / (chi - clo) * 65534) + 1, 0).astype("<u2")
    rgb = np.asarray(Image.open(os.path.join(base, meta["views"]["y"])).convert("RGB")).astype(np.float64) / 255
    lin = rgb ** 2.2
    w = top[..., None]
    col = (lin * w).reshape(m, factor, m, factor, 3).sum((1, 3)) / np.maximum(w.reshape(m, factor, m, factor, 1).sum((1, 3)), 1e-30)
    col = np.clip(col / np.maximum(col.max(-1, keepdims=True), 1e-9), 0, 1)
    col8 = np.round(col * 255).astype(np.uint8)
    return {"m": m, "factor": factor, "range": [clo, chi],
            "data": base64.b64encode(enc.tobytes()).decode(), "colour": base64.b64encode(col8.tobytes()).decode()}


page = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), template)).read()
if "@COARSE@" in page:
    page = page.replace("@COARSE@", json.dumps(coarse()))
for key, value in {"@Y@": data("y"), "@U@": data("x"), "@W@": data("z"), "@STARS@": repr(meta["stars"]),
                   "@FRAME@": json.dumps(meta["frame"]), "@NAME@": meta["name"],
                   "@LY@": light16("y"), "@LU@": light16("x"), "@LW@": light16("z"),
                   "@LIGHTRANGE@": json.dumps(meta.get("light", {}).get("log10Range", [0, 1]))}.items():
    page = page.replace(key, value)
open(out, "w").write(page)
print(f"{out}: {len(page) / 1e6:.1f} MB")
