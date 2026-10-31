"""Topology diagnostics for the BodyParts3D OBJ meshes.

Reports boundary/non-manifold edges before and after welding coincident
vertices, so we know whether the mesh can become a watertight solid.
"""

import os
import sys
from collections import defaultdict


def read_obj(path):
    verts = []
    faces = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("v "):
                p = line.split()
                verts.append((float(p[1]), float(p[2]), float(p[3])))
            elif line.startswith("f "):
                idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[k], idx[k + 1]))
    return verts, faces


def weld(verts, faces, tol):
    keymap = {}
    remap = [0] * len(verts)
    new_verts = []
    q = 1.0 / tol
    for i, (x, y, z) in enumerate(verts):
        key = (int(round(x * q)), int(round(y * q)), int(round(z * q)))
        j = keymap.get(key)
        if j is None:
            j = len(new_verts)
            keymap[key] = j
            new_verts.append((x, y, z))
        remap[i] = j
    new_faces = []
    seen = set()
    for a, b, c in faces:
        t = (remap[a], remap[b], remap[c])
        if len(set(t)) < 3:
            continue
        key = tuple(sorted(t))
        if key in seen:
            continue
        seen.add(key)
        new_faces.append(t)
    return new_verts, new_faces


def edge_report(faces):
    counts = defaultdict(int)
    for a, b, c in faces:
        for e in ((a, b), (b, c), (c, a)):
            counts[(min(e), max(e))] += 1
    boundary = [e for e, n in counts.items() if n == 1]
    nonmanifold = [e for e, n in counts.items() if n > 2]
    return len(counts), boundary, nonmanifold


def boundary_loops(faces):
    _, boundary, _ = edge_report(faces)
    adj = defaultdict(list)
    for a, b in boundary:
        adj[a].append(b)
        adj[b].append(a)
    seen = set()
    loops = 0
    open_chain = 0
    for start in list(adj):
        if start in seen:
            continue
        stack = [start]
        comp = set()
        while stack:
            v = stack.pop()
            if v in comp:
                continue
            comp.add(v)
            stack.extend(adj[v])
        seen |= comp
        degs = [len(adj[v]) for v in comp]
        if all(d == 2 for d in degs):
            loops += 1
        else:
            open_chain += 1
    return loops, open_chain


def main():
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    files = sorted(f for f in os.listdir(folder) if f.lower().endswith(".obj"))
    print(f"{'file':<14}{'V':>7}{'F':>7}{'bnd':>7}{'nonman':>8} | "
          f"{'Vw':>7}{'Fw':>7}{'bnd':>7}{'nonman':>8}{'loops':>7}{'open':>6}  welded?")
    for name in files:
        verts, faces = read_obj(os.path.join(folder, name))
        ne, bnd, nm = edge_report(faces)
        # scale-aware weld tolerance: 1e-4 of the bounding-box diagonal
        xs = [v[0] for v in verts]
        ys = [v[1] for v in verts]
        zs = [v[2] for v in verts]
        diag = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2 +
                (max(zs) - min(zs)) ** 2) ** 0.5
        tol = max(diag * 1e-5, 1e-6)
        wv, wf = weld(verts, faces, tol)
        _, wbnd, wnm = edge_report(wf)
        loops, open_chain = boundary_loops(wf)
        closed = "YES" if (len(wbnd) == 0 and wnm == 0) else "no"
        print(f"{name:<14}{len(verts):>7}{len(faces):>7}{len(bnd):>7}{len(nm):>8} | "
              f"{len(wv):>7}{len(wf):>7}{len(wbnd):>7}{len(wnm):>8}{loops:>7}{open_chain:>6}  {closed}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
