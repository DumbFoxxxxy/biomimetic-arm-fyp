#!/usr/bin/env python3
"""Validate a TPU lattice muscle manufacturing mesh.

Pure numpy (scipy optional, only for speed). No manifold3d / trimesh required.

Checks
  1. connected components of the surface           -> 1 (or extra parts <= --sliver-faces)
  2. closedness of the MAIN component              -> every undirected edge appears twice
  3. wire polyline to material minimum distance    -> must exceed the wire radius
  4. OBJ/STL signed-volume roundtrip               -> relative error < 1e-6

Usage
  python check_tpu_muscle.py <muscle_dir>
  python check_tpu_muscle.py --obj a.obj --stl a.stl --wire design.json
  python check_tpu_muscle.py <dir> --wire-diameter 0.6 --json report.json
  python check_tpu_muscle.py <dir> --sliver-faces 16   # tolerate a tiny debris shell

Exit code 0 = all requested checks passed, 1 = at least one failed.
"""
import argparse
import json
import os
import struct
import sys

import numpy as np


# ---------------------------------------------------------------- loading
def load_obj(path):
    verts, faces = [], []
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('v '):
                p = line.split()
                verts.append((float(p[1]), float(p[2]), float(p[3])))
            elif line.startswith('f '):
                idx = [int(t.split('/')[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[k], idx[k + 1]))
    return np.asarray(verts, dtype=np.float64), np.asarray(faces, dtype=np.int64)


def load_stl(path):
    with open(path, 'rb') as fh:
        head = fh.read(84)
        if len(head) < 84:
            raise ValueError('STL too short: ' + path)
        n = struct.unpack('<I', head[80:84])[0]
        if os.path.getsize(path) != 84 + n * 50:
            raise ValueError('ASCII STL not supported, export binary: ' + path)
        data = np.frombuffer(fh.read(n * 50), dtype=np.uint8).reshape(n, 50)
    tri = data[:, 12:48].copy().view('<f4').reshape(n, 3, 3).astype(np.float64)
    return tri


# ---------------------------------------------------------------- geometry
def weld(verts, faces, tol):
    """Merge coincident vertices. tol <= 0 disables welding (use when indices are already shared)."""
    if tol <= 0:
        return verts, faces
    key = np.round(verts / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return verts[first], inv[faces]


def signed_volume(verts, faces):
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0)


def edge_stats(faces, nverts):
    """Return (edge_id per face-edge matrix, per-edge face counts, boundary, non-manifold)."""
    e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    e.sort(axis=1)
    key = e[:, 0].astype(np.int64) * nverts + e[:, 1]
    _, inverse, counts = np.unique(key, return_inverse=True, return_counts=True)
    return inverse.reshape(3, -1).T, counts


def face_components(edge_ids, nfaces):
    """Union-find over faces sharing an edge; scipy fast path when available."""
    try:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        rows = np.repeat(np.arange(nfaces), 3)
        cols = edge_ids.reshape(-1)
        g = coo_matrix((np.ones(len(rows)), (rows, nfaces + cols)),
                       shape=(nfaces, nfaces + len(np.unique(cols))))
        n, labels = connected_components(g, directed=False)
        return labels
    except Exception:
        pass

    parent = np.arange(nfaces, dtype=np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    flat_edge = edge_ids.reshape(-1)
    flat_face = np.repeat(np.arange(nfaces), 3)
    order = np.argsort(flat_edge, kind='stable')
    fe, ff = flat_edge[order], flat_face[order]
    s = 0
    for i in range(1, len(fe) + 1):
        if i == len(fe) or fe[i] != fe[s]:
            grp = ff[s:i]
            if len(grp) > 1:
                r0 = find(int(grp[0]))
                for g in grp[1:]:
                    r = find(int(g))
                    if r != r0:
                        parent[r] = r0
            s = i
    labels = np.array([find(i) for i in range(nfaces)], dtype=np.int64)
    _, labels = np.unique(labels, return_inverse=True)
    return labels


def point_triangle_distance(p, tri):
    """Closest distance from point p to each triangle in tri (m,3,3)."""
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    ab, ac = b - a, c - a
    ap = p - a
    d1 = np.einsum('ij,ij->i', ab, ap)
    d2 = np.einsum('ij,ij->i', ac, ap)
    bp = p - b
    d3 = np.einsum('ij,ij->i', ab, bp)
    d4 = np.einsum('ij,ij->i', ac, bp)
    cp = p - c
    d5 = np.einsum('ij,ij->i', ab, cp)
    d6 = np.einsum('ij,ij->i', ac, cp)

    closest = np.empty_like(a)
    done = np.zeros(len(tri), dtype=bool)

    m = (d1 <= 0) & (d2 <= 0)
    closest[m] = a[m]; done |= m
    m = (~done) & (d3 >= 0) & (d4 <= d3)
    closest[m] = b[m]; done |= m

    vc = d1 * d4 - d3 * d2
    m = (~done) & (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    den = np.where(np.abs(d1 - d3) < 1e-30, 1e-30, d1 - d3)
    closest[m] = a[m] + (d1[m] / den[m])[:, None] * ab[m]; done |= m

    m = (~done) & (d6 >= 0) & (d5 <= d6)
    closest[m] = c[m]; done |= m

    vb = d5 * d2 - d1 * d6
    m = (~done) & (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    den = np.where(np.abs(d2 - d6) < 1e-30, 1e-30, d2 - d6)
    closest[m] = a[m] + (d2[m] / den[m])[:, None] * ac[m]; done |= m

    va = d3 * d6 - d5 * d4
    m = (~done) & (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    den = np.where(np.abs((d4 - d3) + (d5 - d6)) < 1e-30, 1e-30, (d4 - d3) + (d5 - d6))
    closest[m] = b[m] + ((d4[m] - d3[m]) / den[m])[:, None] * (c[m] - b[m]); done |= m

    m = ~done
    if np.any(m):
        den = va[m] + vb[m] + vc[m]
        den = np.where(np.abs(den) < 1e-30, 1e-30, den)
        closest[m] = a[m] + ab[m] * (vb[m] / den)[:, None] + ac[m] * (vc[m] / den)[:, None]
    return np.linalg.norm(p - closest, axis=1)


def resample(polyline, step):
    pts = np.asarray(polyline, dtype=np.float64)
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        seg = b - a
        n = max(1, int(np.ceil(np.linalg.norm(seg) / step)))
        for k in range(1, n + 1):
            out.append(a + seg * (k / n))
    return np.asarray(out)


def min_clearance(verts, faces, polyline, step, chunk=50000):
    tri = verts[faces]
    pts = resample(polyline, step)
    best = np.inf
    for p in pts:
        for s in range(0, len(tri), chunk):
            d = point_triangle_distance(p, tri[s:s + chunk])
            if d.size:
                best = min(best, float(d.min()))
    return best, len(pts)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('target', nargs='?', help='muscle directory')
    ap.add_argument('--obj')
    ap.add_argument('--stl')
    ap.add_argument('--wire', help='design.json holding wire_polyline_display_mm')
    ap.add_argument('--wire-diameter', type=float, default=0.6)
    ap.add_argument('--step', type=float, default=1.0, help='wire resample step mm')
    ap.add_argument('--weld-tol', type=float, default=0.0,
                    help='vertex weld tolerance mm; 0 = off (manifold3d exports shared indices). '
                         'Enable only when the OBJ duplicates seam vertices.')
    ap.add_argument('--sliver-faces', type=int, default=0,
                    help='tolerate NON-debris extra components up to this face count (default 0)')
    ap.add_argument('--debris-faces', type=int, default=32,
                    help='face count under which an extra part is treated as debris')
    ap.add_argument('--debris-volume', type=float, default=1e-3,
                    help='absolute volume mm3 under which an extra part is treated as debris')
    ap.add_argument('--json', dest='json_out')
    args = ap.parse_args()

    obj, stl, wire = args.obj, args.stl, args.wire
    if args.target:
        base = args.target
        if base.lower().endswith('.obj'):
            obj = obj or base
        elif os.path.isdir(base):
            for attr, name in (('obj', 'TPU_muscle.obj'), ('stl', 'TPU_muscle.stl'), ('wire', 'design.json')):
                cand = os.path.join(base, name)
                if os.path.exists(cand):
                    if attr == 'obj' and not obj:
                        obj = cand
                    elif attr == 'stl' and not stl:
                        stl = cand
                    elif attr == 'wire' and not wire:
                        wire = cand
    if not obj and not stl:
        ap.error('need --obj/--stl or a muscle directory')

    report = {'obj': obj, 'stl': stl, 'wire': wire, 'checks': {}}
    ok = True
    V = F = None

    if obj:
        V, F = load_obj(obj)
        report['checks']['obj_loaded'] = {'verts': int(len(V)), 'faces': int(len(F))}
    if stl:
        report['checks']['stl_loaded'] = {'tris': int(len(load_stl(stl)))}

    if obj:
        Vw, Fw = weld(V, F, args.weld_tol)
        report['checks']['weld'] = {'tolerance_mm': args.weld_tol,
                                    'verts_raw': int(len(V)), 'verts_welded': int(len(Vw))}

        edge_ids, counts = edge_stats(Fw, len(Vw))
        labels = face_components(edge_ids, len(Fw))

        comps = []
        for lab in np.unique(labels):
            sel = labels == lab
            comps.append({'faces': int(sel.sum()),
                          'volume_mm3': round(signed_volume(Vw, Fw[sel]), 6)})
        comps.sort(key=lambda c: c['faces'], reverse=True)
        main_label = np.unique(labels)[np.argmax([c['faces'] for c in comps])]

        def is_debris(c):
            return c['faces'] <= args.debris_faces and abs(c['volume_mm3']) < args.debris_volume

        extras = comps[1:]
        debris = [c for c in extras if is_debris(c)]
        real_extra = [c for c in extras if not is_debris(c)]
        extra_faces = sum(c['faces'] for c in real_extra)

        report['checks']['components'] = {
            'count': len(comps),
            'parts': comps[:8],
            'debris_parts': debris,
            'extra_real_faces': extra_faces,
            'pass': (len(extras) == 0) or (not real_extra) or extra_faces <= args.sliver_faces,
        }
        ok &= (len(extras) == 0) or (not real_extra) or extra_faces <= args.sliver_faces

        main_mask = labels == main_label
        e_main, c_main = edge_stats(Fw[main_mask], len(Vw))
        boundary = int((c_main == 1).sum())
        nonmanifold = int((c_main > 2).sum())
        report['checks']['closed_main'] = {
            'boundary_edges': boundary, 'non_manifold_edges': nonmanifold,
            'pass': boundary == 0 and nonmanifold == 0,
        }
        ok &= (boundary == 0 and nonmanifold == 0)

        all_boundary = int((counts == 1).sum())
        all_nonmanifold = int((counts > 2).sum())
        if all_boundary or all_nonmanifold:
            report['checks']['closed_all'] = {'boundary_edges': all_boundary,
                                              'non_manifold_edges': all_nonmanifold,
                                              'note': 'includes non-main parts'}

        report['checks']['volume_mm3'] = round(signed_volume(V, F), 6)

        if wire:
            with open(wire, 'r', encoding='utf-8') as fh:
                pol = json.load(fh).get('wire_polyline_display_mm')
            if pol:
                dist, npts = min_clearance(V, F, pol, args.step)
                r = args.wire_diameter / 2.0
                report['checks']['wire_clearance'] = {
                    'min_mm': round(dist, 6), 'wire_radius_mm': r,
                    'samples': npts, 'pass': dist > r,
                }
                ok &= dist > r
            else:
                report['checks']['wire_clearance'] = {'skipped': 'no polyline in file'}

    if obj and stl:
        tri = load_stl(stl)
        v_obj = signed_volume(V, F)
        v_stl = signed_volume(tri.reshape(-1, 3), np.arange(tri.shape[0] * 3).reshape(-1, 3))
        err = abs(v_obj - v_stl) / max(abs(v_obj), 1e-12)
        report['checks']['roundtrip'] = {'obj_mm3': v_obj, 'stl_mm3': v_stl,
                                         'relative_error': err, 'pass': err < 1e-6}
        ok &= err < 1e-6

    report['passed'] = bool(ok)
    for name, c in report['checks'].items():
        if isinstance(c, dict) and 'pass' in c:
            print('[%s] %-16s %s' % ('PASS' if c['pass'] else 'FAIL', name, c))
        else:
            print('[INFO] %-16s %s' % (name, c))
    print('\nRESULT:', 'PASS' if ok else 'FAIL')

    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print('report written:', args.json_out)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
