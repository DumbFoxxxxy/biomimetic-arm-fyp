"""Small mesh helpers shared by the tools in this folder (numpy only).

Deliberately minimal: OBJ loading, signed volume, and point/segment distance to a
triangle soup.  The portable muscle validator under .cursor/skills keeps its own
copy so that it stays self-contained.
"""
import numpy as np


def load_obj(path):
    """Return (vertices, triangles) from a Wavefront OBJ."""
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
    return np.asarray(verts, dtype=float), np.asarray(faces, dtype=np.int64)


def signed_volume(verts, faces):
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0)


def point_triangle_distance(p, tri):
    """Distance from point `p` to every triangle in `tri` (m, 3, 3)."""
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


def segment_min_distance(p0, p1, tri, samples=None):
    """Minimum distance from the segment p0-p1 to a triangle soup."""
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    length = float(np.linalg.norm(p1 - p0))
    if samples is None:
        samples = max(2, int(np.ceil(length / 1.0)) + 1)
    best = np.inf
    for t in np.linspace(0.0, 1.0, samples):
        p = p0 + (p1 - p0) * t
        best = min(best, float(point_triangle_distance(p, tri).min()))
    return best
