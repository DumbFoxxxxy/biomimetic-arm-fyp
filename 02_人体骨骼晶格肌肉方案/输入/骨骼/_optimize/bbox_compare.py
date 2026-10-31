"""Compare bounding boxes: raw OBJ vs STL exported from the user's SLDPRT."""

import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decimate import read_obj  # noqa: E402


def read_stl(path):
    with open(path, "rb") as fh:
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        data = fh.read(n * 50)
    arr = np.frombuffer(data, dtype=np.uint8).reshape(n, 50)
    tri = arr[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(np.float64)
    return tri.reshape(-1, 3)


def bbox(V):
    lo = V.min(axis=0)
    hi = V.max(axis=0)
    return lo, hi, hi - lo


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, ".."))
    pairs = [
        ("FJ3349.obj", "\u6861\u9aa8.STL"),
        ("FJ3391.obj", "\u5c3a\u9aa8.STL"),
        ("FJ3368.obj", "\u80b1\u9aa8.STL"),
        ("FJ3384.obj", "\u80a9\u80c4\u9aa8.STL"),
    ]
    for obj, stl in pairs:
        V, _ = read_obj(os.path.join(root, obj))
        S = read_stl(os.path.join(here, "original", stl))
        lo1, hi1, d1 = bbox(V)
        lo2, hi2, d2 = bbox(S)
        print(f"--- {obj} vs {stl}")
        print(f"  obj lo={np.round(lo1,2)} hi={np.round(hi1,2)} dim={np.round(d1,2)}")
        print(f"  stl lo={np.round(lo2,2)} hi={np.round(hi2,2)} dim={np.round(d2,2)}")
        print(f"  dim ratio obj/stl = {np.round(d1 / d2, 4)}")


if __name__ == "__main__":
    main()
