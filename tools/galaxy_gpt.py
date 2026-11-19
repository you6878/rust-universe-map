#!/usr/bin/env python3
"""A galaxy's three pictures for galaxy_xyz.py, made by OpenAI's image model (gpt-image-2): the view
along its shortest axis (from above, Y) and the two side views along the other two (X, Z), the side
views drawn from the top one so that they are the same galaxy.

Sizes (memory: galaxy-picture-sizes): the top view at 2048 (the near picture, twice the X·Y·Z size),
the side views at 1024. The X·Y·Z 1024 and the far 512 are made from the top one later.

Safety:
  - the key is read from the environment (OPENAI_API_KEY) or --key-file (a line KEY=VALUE or the bare
    key); it is never printed, logged or put on a command line;
  - --dry-run prints what would be asked (prompts, sizes, number of calls) and asks nothing;
  - a picture already in OUT is not asked for again;
  - every reply's usage is written to OUT/usage.json.

Usage: galaxy_gpt.py OUT_DIR --kind elliptical|spiral|irregular [--name N] [--axis-ratio 0.8]
                     [--quality low|medium|high] [--key-file F] [--key-var GPT_IMAGE_API_KEY] [--dry-run]
"""
import argparse
import base64
import json
import os
import sys
import time
import urllib.request
import uuid

API = "https://api.openai.com/v1/images"
MODEL = "gpt-image-2"

args = argparse.ArgumentParser()
args.add_argument("out")
args.add_argument("--kind", choices=["elliptical", "spiral", "irregular"], required=True)
args.add_argument("--name", default="")
args.add_argument("--axis-ratio", type=float, default=0.8, help="seen from above: short over long")
args.add_argument("--thickness", type=float, default=0.6, help="seen from the side: height over width")
args.add_argument("--quality", choices=["low", "medium", "high"], default="low")
args.add_argument("--key-file")
args.add_argument("--key-var", default="OPENAI_API_KEY")
args.add_argument("--dry-run", action="store_true")
a = args.parse_args()

TOP = {
    "elliptical": "a giant elliptical galaxy at the centre of a galaxy cluster, seen straight along its shortest axis: a smooth, "
                  "featureless warm golden-white glow, axis ratio about {q:.1f}, a bright compact core fading smoothly outward "
                  "(de Vaucouleurs profile), a faint extended halo, a few faint globular clusters as tiny points near it",
    "spiral": "a spiral galaxy seen exactly face-on: a yellowish bulge, blue spiral arms with pink star-forming regions and "
              "brown dust lanes along their inner edges",
    "irregular": "an irregular galaxy seen from above: clumpy, blue, knots of young stars and pink gas, no clear arms",
}
SIDE = {
    "elliptical": "a triaxial elliptical, not a disc: no dust lane, no arms; height over width about {t:.1f}",
    "spiral": "a thin disc seen exactly edge-on with a dust lane along its middle and the bulge rising above and below it",
    "irregular": "a thick clumpy body seen from the side, height over width about {t:.1f}",
}
COMMON = "Astronomical photograph, exactly centred, the galaxy filling about 80% of the frame, pure black background, " \
         "no foreground stars, no other galaxies, no text, no labels."

top_prompt = f"{TOP[a.kind].format(q=a.axis_ratio)[0].upper()}{TOP[a.kind].format(q=a.axis_ratio)[1:]}. {COMMON}"
sides = {
    "x": "its longest axis lies left to right",
    "z": "its intermediate axis lies left to right",
}
side_prompts = {v: (f"This is a galaxy seen from above. Draw the SAME galaxy seen from the side, exactly perpendicular to the "
                    f"first view ({how}): the same colours, the same core and light profile, the same extent left to right, as "
                    f"{SIDE[a.kind].format(t=a.thickness)}. The same scale as the input. {COMMON}") for v, how in sides.items()}

plan = [("y", "generations", "2048x2048", top_prompt)] + [(v, "edits", "1024x1024", p) for v, p in side_prompts.items()]
todo = [p for p in plan if not os.path.exists(os.path.join(a.out, f"{p[0]}.png"))]
print(f"{a.name or a.kind}: {len(todo)} of {len(plan)} pictures to ask {MODEL} for, quality {a.quality}")
for v, kind, size, prompt in plan:
    there = "" if any(t[0] == v for t in todo) else " (there already)"
    print(f"  {v} {kind} {size}{there}: {prompt}")
if a.dry_run or not todo:
    sys.exit(0)


def read_key() -> str:
    if a.key_file:
        for line in open(a.key_file):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k, eq, v = line.partition("=")
            if not eq:
                return k.strip().strip("'\"")
            if k.strip() == a.key_var:
                return v.strip().strip("'\"")
        sys.exit(f"{a.key_file}: no {a.key_var}")
    key = os.environ.get(a.key_var, "")
    if not key:
        sys.exit(f"no key: set {a.key_var} or give --key-file")
    return key


KEY = read_key()


def ask(kind: str, body: bytes, content_type: str) -> dict:
    req = urllib.request.Request(f"{API}/{kind}", data=body, method="POST",
                                 headers={"Authorization": f"Bearer {KEY}", "Content-Type": content_type})
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        try:
            msg = json.load(e).get("error", {}).get("message", "")
        except Exception:
            msg = ""
        sys.exit(f"{kind}: HTTP {e.code}: {msg[:300]}")


def multipart(fields: dict, files: dict) -> tuple:
    b = uuid.uuid4().hex
    out = bytearray()
    for k, v in fields.items():
        out += f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode()
    for k, (fname, data) in files.items():
        out += f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"; filename=\"{fname}\"\r\nContent-Type: image/png\r\n\r\n".encode()
        out += data + b"\r\n"
    out += f"--{b}--\r\n".encode()
    return bytes(out), f"multipart/form-data; boundary={b}"


os.makedirs(a.out, exist_ok=True)
usage_file = os.path.join(a.out, "usage.json")
usage = json.load(open(usage_file)) if os.path.exists(usage_file) else []
for v, kind, size, prompt in todo:
    t = time.time()
    if kind == "generations":
        reply = ask(kind, json.dumps({"model": MODEL, "prompt": prompt, "size": size, "quality": a.quality, "n": 1}).encode(), "application/json")
    else:
        top = open(os.path.join(a.out, "y.png"), "rb").read()
        body, ctype = multipart({"model": MODEL, "prompt": prompt, "size": size, "quality": a.quality}, {"image[]": ("y.png", top)})
        reply = ask(kind, body, ctype)
    open(os.path.join(a.out, f"{v}.png"), "wb").write(base64.b64decode(reply["data"][0]["b64_json"]))
    usage.append({"view": v, "size": size, "quality": a.quality, "seconds": round(time.time() - t, 1), "usage": reply.get("usage")})
    json.dump(usage, open(usage_file, "w"), indent=1)
    print(f"  {v}: {size} in {time.time() - t:.0f} s, {reply.get('usage', {}).get('output_tokens')} output tokens")
