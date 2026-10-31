"""Quality checks that SolidWorks cares about before knitting a solid."""

import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decimate import read_obj, weld  # noqa: E402

import numpy as np  # noqa: E402


def check(V, F, label):
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    area = 0.5 * np.linalg.norm(n, axis=1)

    degenerate = int((area <= 1e-12).sum())
    tiny = int((area <= 1e-8).sum())

    seen = defaultdict(int)
    for t in F:
        seen[tuple(sorted(int(x) for x in t))] += 1
    dup = sum(1 for v in seen.values() if v > 1)

    # orientation: interior edges must be traversed in opposite directions
    dir_count = defaultdict(int)
    for a, b, c in F:
        dir_count[(int(a), int(b))] += 1
        dir_count[(int(b), int(c))] += 1
        dir_count[(int(c), int(a))] += 1
    bad_orient = 0
    for (a, b), cnt in dir_count.items():
        if cnt > 1:
            bad_orient += cnt - 1

    d = np.linalg.norm(V.max(axis=0) - V.min(axis=0))
    print(f"{label:<12} V={len(V):>6} F={len(F):>6} degen={degenerate:>5} tiny={tiny:>5} "
          f"dupFaces={dup:>4} badOrientEdges={bad_orient:>5} diag={d:.2f} "
          f"minArea={area.min():.3e}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, ".."))
    for name in sorted(os.listdir(root)):
        if not name.lower().endswith(".obj"):
            continue
        V, F = read_obj(os.path.join(root, name))
        diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
        Vw, Fw = weld(V, F, max(diag * 1e-5, 1e-6))
        check(V, F, name + ":raw")
        check(Vw, Fw, name + ":weld")


if __name__ == "__main__":
    main()
