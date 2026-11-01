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


def _closest_all(p, tri):
    """Closest point on every triangle to `p`: (distances, closest points)."""
    p = np.asarray(p, float)
    tri = np.asarray(tri, float)
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

    return np.linalg.norm(p - closest, axis=1), closest


def point_triangle_distance(p, tri):
    """Distance from point `p` to every triangle in `tri` (m, 3, 3)."""
    return _closest_all(np.asarray(p, float), np.asarray(tri, float))[0]


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


def signed_distance(p, tri):
    """Signed distance from `p` to a closed triangle mesh.

    Negative inside, positive outside.  The sign comes from the nearest
    triangle's own normal, so it is reliable except within one triangle of a
    crease or a very thin wall - which is exactly the regime that has to be
    reported rather than trusted.
    """
    p = np.asarray(p, float)
    tri = np.asarray(tri, float)
    dist, closest = _closest_all(p, tri)
    k = int(np.argmin(dist))
    d = float(dist[k])
    n = np.cross(tri[k, 1] - tri[k, 0], tri[k, 2] - tri[k, 0])
    norm = np.linalg.norm(n)
    if norm < 1e-12 or d < 1e-12:
        return d
    n = n / norm
    return -d if (p - closest[k]) @ n < 0.0 else d


def edge_census(faces):
    """(edges, boundary edges, non-manifold edges) of a triangle list."""
    f = np.asarray(faces, dtype=np.int64)
    e = np.sort(np.c_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]].reshape(-1, 2), axis=1)
    uniq, counts = np.unique(e, axis=0, return_counts=True)
    return len(uniq), int((counts == 1).sum()), int((counts > 2).sum())


def weld(verts, faces, decimals=4):
    """Merge coincident vertices and drop degenerate triangles.

    The BodyParts3D OBJ exports carry duplicate vertices per face, so the raw
    soup reports hundreds of boundary edges and ``manifold3d`` silently returns
    an EMPTY solid for it - every interference query against such a bone then
    answers 0.  Welding by rounded position closes all three arm bones exactly
    (humerus / scapula / clavicle: 0 boundary edges, 0 non-manifold edges), so
    always weld before building a bone solid.

    Check with ``edge_census``: the welded result must report 0 boundary edges.
    """
    V = np.asarray(verts, dtype=float)
    F = np.asarray(faces, dtype=np.int64)
    if len(V) == 0 or len(F) == 0:
        return V, F
    _, first, inverse = np.unique(np.round(V, decimals), axis=0,
                                 return_index=True, return_inverse=True)
    newF = inverse.reshape(-1)[F.reshape(-1)].reshape(F.shape)
    keep = ((newF[:, 0] != newF[:, 1]) & (newF[:, 1] != newF[:, 2])
            & (newF[:, 0] != newF[:, 2]))
    return V[first], newF[keep]


# ------------------------------------------------------------------ self test
def _split_cube(size=10.0):
    """A cube as a triangle soup with a *separate* vertex copy per triangle."""
    h = size / 2.0
    corners = np.array([[-h, -h, -h], [h, -h, -h], [h, h, -h], [-h, h, -h],
                        [-h, -h, h], [h, -h, h], [h, h, h], [-h, h, h]])
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (2, 3, 7, 6), (1, 2, 6, 5), (0, 4, 7, 3)]
    V, F = [], []
    for a, b, c, d in quads:
        for tri in ((a, b, c), (a, c, d)):
            base = len(V)
            V.extend(corners[list(tri)])
            F.append((base, base + 1, base + 2))
    return np.asarray(V, float), np.asarray(F, np.int64)


def self_test():
    ok = True

    def check(name, value, expected, tol=1e-9):
        nonlocal ok
        good = abs(value - expected) <= tol
        print('[%s] %-38s %.6f (expect %.6f)'
              % ('PASS' if good else 'FAIL', name, value, expected))
        ok &= good

    V, F = _split_cube(10.0)
    e, b, n = edge_census(F)
    check('split cube boundary edges', b, 36)
    W, G = weld(V, F)
    check('welded vertex count', len(W), 8)
    check('welded triangle count', len(G), 12)
    e2, b2, n2 = edge_census(G)
    check('welded boundary edges', b2, 0)
    check('welded non-manifold edges', n2, 0)
    check('welded volume', abs(signed_volume(W, G)), 1000.0)
    print('\nRESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    import sys
    sys.exit(self_test())
