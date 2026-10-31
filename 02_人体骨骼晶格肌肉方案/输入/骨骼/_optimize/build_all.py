"""Weld + decimate the four mesh bones exported from the user's SLDPRT files.

Source : ./original/<bone>.STL   (exact coordinates the user already has)
Output : ./reduced/<ascii>.stl   (clean, watertight, ready to import as a solid)
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decimate import read_stl, weld, decimate, write_stl  # noqa: E402

# (source name, ascii output name, target triangle count)
JOBS = [
    ("\u6861\u9aa8.STL", "radius", 600),     # 桡骨 radius
    ("\u5c3a\u9aa8.STL", "ulna", 700),       # 尺骨 ulna
    ("\u80b1\u9aa8.STL", "humerus", 1200),   # 肱骨 humerus
    ("\u80a9\u80db\u9aa8.STL", "scapula", 2500),  # 肩胛骨 scapula
]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.join(here, "original")
    out_dir = os.path.join(here, "reduced")
    os.makedirs(out_dir, exist_ok=True)
    for src_name, out_name, target in JOBS:
        src = os.path.join(src_dir, src_name)
        if not os.path.exists(src):
            print(f"missing {src}")
            continue
        t0 = time.time()
        V, F = read_stl(src)
        diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
        Vw, Fw = weld(V, F, max(diag * 1e-6, 1e-7))
        print(f"{out_name}: stl-weld V={len(Vw)} F={len(Fw)}  ({time.time()-t0:.1f}s)")
        sys.stdout.flush()
        t1 = time.time()
        Vr, Fr = decimate(Vw, Fw, target)
        out = os.path.join(out_dir, out_name + ".stl")
        write_stl(out, Vr, Fr)
        print(f"  decimated -> F={len(Fr)} V={len(Vr)}  ({time.time()-t1:.1f}s)  "
              f"bbox={np.round(Vr.min(axis=0),2)}..{np.round(Vr.max(axis=0),2)}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
