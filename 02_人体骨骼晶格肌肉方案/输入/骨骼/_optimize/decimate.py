"""Weld + quadric-error decimation for the BodyParts3D bone meshes.

Input : .obj meshes in the parent folder
Output: clean, watertight binary .stl files in ./welded and ./reduced

The OBJ files ship with unwelded seams; welding first is what turns them
into closed manifold shells that SolidWorks can knit into a solid body.
"""

import argparse
import heapq
import os
import struct

import numpy as np


def read_stl(path):
    with open(path, "rb") as fh:
        head = fh.read(84)
    if head[:5].lower() == b"solid" and b"facet" in head:
        return read_stl_ascii(path)
    with open(path, "rb") as fh:
        fh.seek(80)
        n = struct.unpack("<I", fh.read(4))[0]
        data = fh.read(n * 50)
    arr = np.frombuffer(data, dtype=np.uint8).reshape(n, 50)
    tri = arr[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(np.float64)
    return weld_stl(tri)


def read_stl_ascii(path):
    tris = []
    cur = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("vertex"):
                p = s.split()
                cur.append((float(p[1]), float(p[2]), float(p[3])))
                if len(cur) == 3:
                    tris.append(tuple(cur))
                    cur = []
    return weld_stl(np.asarray(tris, dtype=np.float64))


def weld_stl(tris):
    """STL has no shared indices; rebuild them by welding coincident corners."""
    flat = tris.reshape(-1, 3)
    keymap = {}
    remap = np.zeros(len(flat), dtype=np.int64)
    verts = []
    for i, p in enumerate(flat):
        key = (repr(round(p[0], 6)), repr(round(p[1], 6)), repr(round(p[2], 6)))
        j = keymap.get(key)
        if j is None:
            j = len(verts)
            keymap[key] = j
            verts.append(p)
        remap[i] = j
    faces = remap.reshape(-1, 3)
    keep = []
    seen = set()
    for t in faces:
        if len(set(int(x) for x in t)) < 3:
            continue
        key = tuple(sorted(int(x) for x in t))
        if key in seen:
            continue
        seen.add(key)
        keep.append(t)
    return np.asarray(verts, dtype=np.float64), np.asarray(keep, dtype=np.int64)


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
    return np.asarray(verts, dtype=np.float64), np.asarray(faces, dtype=np.int64)


def read_mesh(path):
    if path.lower().endswith(".stl"):
        return read_stl(path)
    return read_obj(path)


def weld(verts, faces, tol):
    q = 1.0 / tol
    keymap = {}
    remap = np.zeros(len(verts), dtype=np.int64)
    new_verts = []
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
        t = (int(remap[a]), int(remap[b]), int(remap[c]))
        if len(set(t)) < 3:
            continue
        key = tuple(sorted(t))
        if key in seen:
            continue
        seen.add(key)
        new_faces.append(t)
    return np.asarray(new_verts, dtype=np.float64), np.asarray(new_faces, dtype=np.int64)


def face_planes(V, F):
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    n = np.cross(p1 - p0, p2 - p0)
    norm = np.linalg.norm(n, axis=1)
    area = 0.5 * norm
    ok = norm > 1e-14
    nrm = np.zeros_like(n)
    nrm[ok] = n[ok] / norm[ok, None]
    d = -np.einsum("ij,ij->i", nrm, p0)
    return nrm, d, area


def edge_map(F):
    counts = {}
    for a, b, c in F:
        for e in ((a, b), (b, c), (c, a)):
            k = (min(e), max(e))
            counts[k] = counts.get(k, 0) + 1
    return counts


def build_quadrics(V, F, boundary_weight):
    nrm, d, area = face_planes(V, F)
    Q = np.zeros((len(V), 4, 4))
    K = np.zeros((len(F), 4, 4))
    K[:, :3, :3] = nrm[:, :, None] * nrm[:, None, :]
    K[:, :3, 3] = nrm * d[:, None]
    K[:, 3, :3] = nrm * d[:, None]
    K[:, 3, 3] = d * d
    K *= area[:, None, None]
    for j in range(3):
        np.add.at(Q, F[:, j], K)

    counts = edge_map(F)
    bedges = [k for k, v in counts.items() if v == 1]
    if bedges:
        edge_face = {}
        for fi in range(len(F)):
            t = F[fi]
            for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                k = (min(int(e[0]), int(e[1])), max(int(e[0]), int(e[1])))
                edge_face.setdefault(k, fi)
        kf = np.zeros((len(bedges), 4, 4))
        for i, (a, b) in enumerate(bedges):
            e = V[b] - V[a]
            el = np.linalg.norm(e)
            if el < 1e-14:
                continue
            e = e / el
            # use any incident face normal to build the boundary plane
            fi = edge_face.get((min(a, b), max(a, b)))
            fn = None
            if fi is not None:
                t = F[fi]
                fn = np.cross(V[t[1]] - V[t[0]], V[t[2]] - V[t[0]])
            if fn is None or np.linalg.norm(fn) < 1e-14:
                continue
            bn = np.cross(fn / np.linalg.norm(fn), e)
            nb = np.linalg.norm(bn)
            if nb < 1e-14:
                continue
            bn = bn / nb
            dd = -bn.dot(V[a])
            kf[i, :3, :3] = boundary_weight * np.outer(bn, bn)
            kf[i, :3, 3] = boundary_weight * bn * dd
            kf[i, 3, :3] = boundary_weight * bn * dd
            kf[i, 3, 3] = boundary_weight * dd * dd
        for i, (a, b) in enumerate(bedges):
            Q[a] += kf[i]
            Q[b] += kf[i]
    return Q


def decimate(V, F, target_faces, boundary_weight=1000.0, max_edge_factor=10.0):
    V = np.array(V, dtype=np.float64, copy=True)
    F = np.array(F, dtype=np.int64, copy=True)
    nv = len(V)
    Q = build_quadrics(V, F, boundary_weight)

    vert_faces = [set() for _ in range(nv)]
    for fi in range(len(F)):
        for j in range(3):
            vert_faces[F[fi, j]].add(fi)

    face_alive = np.ones(len(F), dtype=bool)
    vert_alive = np.ones(nv, dtype=bool)
    version = np.zeros(nv, dtype=np.int64)

    def solve_edge(u, v):
        Qb = Q[u] + Q[v]
        A = Qb[:3, :3]
        rhs = -Qb[:3, 3]
        mid = 0.5 * (V[u] + V[v])
        x = None
        try:
            if abs(np.linalg.det(A)) > 1e-18:
                cand = np.linalg.solve(A, rhs)
                if np.all(np.isfinite(cand)):
                    x = cand
        except np.linalg.LinAlgError:
            x = None
        if x is None or np.linalg.norm(x - mid) > max_edge_factor * np.linalg.norm(V[v] - V[u]) + 1e-12:
            x = mid
        w = np.array([x[0], x[1], x[2], 1.0])
        return float(w @ Qb @ w), x

    counts = edge_map(F)
    heap = []
    for (a, b) in counts:
        c, x = solve_edge(a, b)
        heap.append((c, int(a), int(b), 0, 0, x[0], x[1], x[2]))
    heapq.heapify(heap)

    alive_faces = int(face_alive.sum())
    while heap and alive_faces > target_faces:
        c, u, v, vu, vv, x0, x1, x2 = heapq.heappop(heap)
        if not (vert_alive[u] and vert_alive[v]):
            continue
        if version[u] != vu or version[v] != vv:
            continue
        common = vert_faces[u] & vert_faces[v]
        if len(common) != 2:
            continue  # link condition: keep the surface manifold

        newpos = np.array([x0, x1, x2])
        # reject collapses that would flip adjacent face normals
        flip = False
        for fi in vert_faces[u] | vert_faces[v]:
            if not face_alive[fi]:
                continue
            t = F[fi]
            p = [newpos if w in (u, v) else V[w] for w in t]
            na = np.cross(V[t[1]] - V[t[0]], V[t[2]] - V[t[0]])
            nb = np.cross(p[1] - p[0], p[2] - p[0])
            la, lb = np.linalg.norm(na), np.linalg.norm(nb)
            if la > 1e-14 and lb > 1e-14 and na.dot(nb) / (la * lb) < 0.0:
                flip = True
                break
        if flip:
            continue

        V[u] = newpos
        Q[u] = Q[u] + Q[v]
        for fi in list(vert_faces[v]):
            if not face_alive[fi]:
                continue
            row = F[fi]
            row[row == v] = u
            if row[0] == row[1] or row[1] == row[2] or row[0] == row[2]:
                face_alive[fi] = False
                alive_faces -= 1
                for w in set(int(t) for t in row):
                    vert_faces[w].discard(fi)
            else:
                vert_faces[u].add(fi)
        vert_faces[v] = set()
        vert_alive[v] = False
        version[u] += 1
        version[v] += 1

        for fi in list(vert_faces[u]):
            if not face_alive[fi]:
                continue
            for w in set(int(t) for t in F[fi]):
                if w != u and vert_alive[w]:
                    cc, xx = solve_edge(u, w)
                    heapq.heappush(heap, (cc, u, w, version[u], version[w],
                                          xx[0], xx[1], xx[2]))

    F2 = F[face_alive]
    used = np.unique(F2)
    remap = -np.ones(nv, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return V[used], remap[F2]


def write_stl(path, V, F):
    with open(path, "wb") as fh:
        fh.write(b"\0" * 80)
        fh.write(struct.pack("<I", len(F)))
        p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
        n = np.cross(p1 - p0, p2 - p0)
        ln = np.linalg.norm(n, axis=1)
        ok = ln > 1e-14
        n[ok] /= ln[ok, None]
        n[~ok] = 0.0
        rec = np.zeros((len(F), 12), dtype=np.float32)
        rec[:, 0:3] = n
        rec[:, 3:6] = p0
        rec[:, 6:9] = p1
        rec[:, 9:12] = p2
        rec = rec.astype("<f4")
        for r in rec:
            fh.write(r.tobytes())
            fh.write(b"\0\0")


def write_obj(path, V, F):
    with open(path, "w") as fh:
        for v in V:
            fh.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for f in F:
            fh.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=0,
                    help="target triangle count; 0 = weld only")
    ap.add_argument("--src", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(here, ".."))
    src = args.src or os.path.join(root, "FJ3384.obj")
    V, F = read_mesh(src)
    diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
    Vw, Fw = weld(V, F, max(diag * 1e-5, 1e-6))
    print(f"weld: V {len(V)} -> {len(Vw)},  F {len(F)} -> {len(Fw)}")
    Vr, Fr = Vw, Fw
    if args.target and args.target < len(Fw):
        Vr, Fr = decimate(Vw, Fw, args.target)
        print(f"decimate: F {len(Fw)} -> {len(Fr)},  V {len(Vw)} -> {len(Vr)}")
    out = args.out or os.path.join(here, "out.stl")
    write_stl(out, Vr, Fr)
    print(f"wrote {out}  tris={len(Fr)}")


if __name__ == "__main__":
    main()
