#!/usr/bin/env python3
"""Measured heights of a planet, as the grid of numbers the viewer's planet object reads.

    make_terrain.py SOURCE OUT_DIR NAME [--width 10800] [--scale-km KM] [--preview]

SOURCE is either
  * a NetCDF file of elevations in metres on a grid of longitude and latitude, such as NOAA's
    ETOPO 2022 (variable `z`; land heights and sea depths), or
  * a picture whose brightness is the height (the height maps of the Gaia Sky data); then
    --scale-km says how high white is. Such a picture has no sea floor: the sea is at 0.

Written to OUT_DIR:
  NAME.i16.gz        heights in whole metres, 16 bits little-endian, row by row from the north
                     pole to the south, each row eastwards from longitude -180
  NAME_metadata.txt  "Grid size: W H" and the lowest and highest height

The grid is made smaller than the source by taking the mean of the source points in each cell.
"""
import argparse
import gzip
import os
import sys

import numpy as np


def read_netcdf(path):
    import h5py
    with h5py.File(path, "r") as f:
        z = f["z"]
        lat = f["lat"][:]
        heights = np.empty(z.shape, dtype=np.float32)
        z.read_direct(heights)
    if lat[0] < lat[-1]:          # south first: turn it over
        heights = heights[::-1]
    return heights


def read_picture(path, scale_km):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    a = np.asarray(Image.open(path).convert("L"), dtype=np.float32)
    return a / 255.0 * scale_km * 1000.0


def shrink(heights, width):
    h, w = heights.shape
    if width >= w:
        return heights
    if w % width == 0 and h % (width // 2) == 0:
        k = w // width
        return heights.reshape(h // k, k, width, k).mean(axis=(1, 3))
    from scipy.ndimage import zoom
    return zoom(heights, (width / 2 / h, width / w), order=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("out_dir")
    ap.add_argument("name")
    ap.add_argument("--width", type=int, default=10800)
    ap.add_argument("--scale-km", type=float, default=None)
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()
    if args.source.lower().endswith((".nc", ".nc4")):
        heights = read_netcdf(args.source)
    else:
        if args.scale_km is None:
            sys.exit("a picture needs --scale-km: how high white is")
        heights = read_picture(args.source, args.scale_km)
    heights = shrink(heights, args.width)
    grid = np.clip(np.rint(heights), -32768, 32767).astype("<i2")
    h, w = grid.shape
    os.makedirs(args.out_dir, exist_ok=True)
    out = os.path.join(args.out_dir, args.name)
    with gzip.open(out + ".i16.gz", "wb", compresslevel=6) as f:
        f.write(grid.tobytes())
    with open(out + "_metadata.txt", "w") as f:
        f.write(f"Grid size: {w} {h}\nHeights in metres: {grid.min()} {grid.max()}\nSource: {os.path.basename(args.source)}\n")
    print(f"{out}.i16.gz: {w} x {h}, {grid.min()} to {grid.max()} m, {os.path.getsize(out + '.i16.gz') / 1e6:.1f} MB")
    if args.preview:
        from PIL import Image
        small = grid[:: max(h // 900, 1), :: max(w // 1800, 1)].astype(np.float32)
        land = np.clip(small / 6000.0, 0, 1)[..., None] * np.array([200, 160, 120]) + np.array([40, 90, 40])
        sea = np.clip(1 + small / 8000.0, 0, 1)[..., None] * np.array([40, 110, 200])
        Image.fromarray(np.where(small[..., None] > 0, land, sea).clip(0, 255).astype(np.uint8)).save(out + "_preview.jpg", quality=85)


if __name__ == "__main__":
    main()
