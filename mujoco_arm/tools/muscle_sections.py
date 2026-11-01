"""Planar sections through a muscle solid, and 2D metrics on them.

Two studies need the same thing - the muscle envelope study (does a channel fit?)
and the end-cap study (does a rigid cap fit?) - so the slicing lives here.

A section is taken perpendicular to the origin->insertion chord, in a frame whose
+z is that chord.  The chord is only a *reference direction*, not a fibre
direction and not a tendon path: both of its ends sit on the bone surface, which
is a boundary of the belly, so the chord leaves the muscle almost immediately.
That is exactly why the interesting quantity is the section centroid curve and the
wall from it to the section boundary, not the distance along the chord.

All distances are mm, everything stays in source coordinates unless a name ends
in ``_xy`` (that is the slicing frame, 2D).
"""
import numpy as np

import frames

#: grid step used for the largest-inscribed-circle search [mm]
INRADIUS_STEP_MM = 1.0


def rot_to_z(direction):
    """Rotation matrix whose third column is `direction` (chord -> +z)."""
    z = frames.unit(direction)
    helper = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    x = frames.unit(np.cross(helper, z))
    y = np.cross(z, x)
    return np.stack([x, y, z], axis=1)


def _edges(polys):
    """All polygon edges as (m, 2, 2) arrays."""
    out = []
    for poly in polys:
        p = np.asarray(poly, float)
        if len(p) < 2:
            continue
        out.append(np.stack([p, np.roll(p, -1, axis=0)], axis=1))
    if not out:
        return np.zeros((0, 2, 2))
    return np.concatenate(out, axis=0)


def _segment_distance(points, edges):
    """Distance from every point to the nearest edge, vectorised per edge."""
    best = np.full(len(points), np.inf)
    for a, b in edges:
        ab = b - a
        denom = float(ab @ ab)
        if denom < 1e-18:
            d = np.linalg.norm(points - a, axis=1)
        else:
            t = np.clip(((points - a) @ ab) / denom, 0.0, 1.0)
            d = np.linalg.norm(points - (a + t[:, None] * ab), axis=1)
        best = np.minimum(best, d)
    return best


def _inside_even_odd(points, edges):
    """Even-odd rule over all edges: holes and multiple loops both work."""
    inside = np.zeros(len(points), dtype=bool)
    for a, b in edges:
        ay, by = a[1], b[1]
        if ay == by:
            continue
        crosses = (ay > 0.0) != (by > 0.0)
        if not np.any(crosses):
            continue
        dy = by - ay
        x_at = a[0] + (b[0] - a[0]) * (-ay) / dy
        inside ^= crosses & (points[:, 0] < x_at)
    return inside


def _area_centroid(polys):
    """Signed-area weighted centroid; holes enter with the opposite sign."""
    total = 0.0
    acc = np.zeros(2)
    for poly in polys:
        p = np.asarray(poly, float)
        if len(p) < 3:
            continue
        x, y = p[:, 0], p[:, 1]
        x1, y1 = np.roll(x, -1), np.roll(y, -1)
        cross = x * y1 - x1 * y
        a = 0.5 * cross.sum()
        if abs(a) < 1e-12:
            continue
        total += a
        acc += a * np.array([((x + x1) * cross).sum(), ((y + y1) * cross).sum()]) / (6.0 * a)
    if abs(total) < 1e-12:
        return 0.0, None
    return total, acc / total


def path_clearance(points, tri, step=1.5, wall_mm=1.6):
    """Wall along a polyline inside a closed mesh.

    Sign convention, which is easy to get backwards: ``mesh_geom.signed_distance``
    returns NEGATIVE inside the material.  So along a route

        worst_mm    = max(signed)  the least-inside point  (the binding one)
        deepest_mm  = min(signed)  the most-inside point

    A route can carry a channel only when even its worst point is at least `wall_mm`
    inside, i.e. ``worst_mm <= -wall_mm``.  A positive ``worst_mm`` means the route
    leaves the belly.
    """
    import mesh_geom
    samples = []
    for a, b in zip(points, points[1:]):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        n = max(2, int(np.ceil(np.linalg.norm(b - a) / step)) + 1)
        for t in np.linspace(0.0, 1.0, n):
            samples.append(a + (b - a) * t)
    d = [mesh_geom.signed_distance(p, tri) for p in samples]
    worst = float(max(d))
    return {'samples': len(d),
            'worst_mm': round(worst, 3),
            'deepest_mm': round(float(min(d)), 3),
            'inside_fraction': round(sum(1 for x in d if x < 0) / max(1, len(d)), 3),
            'required_wall_mm': wall_mm,
            'fits_channel': bool(worst <= -wall_mm),
            'profile_mm': [round(float(x), 2) for x in d]}


def section_at(solid, origin, insertion, t, step=INRADIUS_STEP_MM):
    """One perpendicular section at fraction `t` of the origin->insertion chord.

    Returns area, centroid (2D and in source mm), the wall from the centroid to
    the boundary, and the largest inscribed circle found on a `step` grid.
    `inradius_mm` is a lower bound: a 1 mm grid under-reports by up to ~0.7 mm.
    """
    origin = np.asarray(origin, float)
    insertion = np.asarray(insertion, float)
    chord = insertion - origin
    length = float(np.linalg.norm(chord))
    R = rot_to_z(chord)
    moved = solid.transform(np.c_[R.T, -R.T @ origin])
    cs = moved.slice(float(t * length))
    if cs.is_empty():
        return None
    polys = cs.to_polygons()
    area, cen = _area_centroid(polys)
    if cen is None:
        return {'t': t, 'z_mm': t * length, 'area_mm2': abs(cs.area()),
                'centroid_xy_mm': None, 'chord_offset_mm': None,
                'wall_mm': None, 'inradius_mm': None, 'polys_xy': polys}
    edges = _edges(polys)
    wall = float(_segment_distance(cen[None, :], edges)[0])
    lo = edges.reshape(-1, 2).min(axis=0)
    hi = edges.reshape(-1, 2).max(axis=0)
    nx = max(1, int(np.ceil((hi[0] - lo[0]) / step)) + 1)
    ny = max(1, int(np.ceil((hi[1] - lo[1]) / step)) + 1)
    gx, gy = np.meshgrid(lo[0] + step * np.arange(nx), lo[1] + step * np.arange(ny))
    grid = np.c_[gx.ravel(), gy.ravel()]
    inside = _inside_even_odd(grid, edges)
    inradius = None
    inradius_xy = None
    if np.any(inside):
        d = _segment_distance(grid[inside], edges)
        k = int(np.argmax(d))
        inradius = float(d[k])
        inradius_xy = grid[inside][k]
    world = origin + R @ np.array([cen[0], cen[1], t * length])
    world_in = (None if inradius_xy is None else
                origin + R @ np.array([inradius_xy[0], inradius_xy[1], t * length]))
    return {'t': float(t), 'z_mm': round(t * length, 2),
            'area_mm2': round(abs(area), 2),
            'centroid_xy_mm': np.round(cen, 3).tolist(),
            'chord_offset_mm': round(float(np.linalg.norm(cen)), 3),
            'wall_mm': round(wall, 3),
            'inradius_mm': None if inradius is None else round(inradius, 3),
            'inradius_xy_mm': None if inradius_xy is None else np.round(inradius_xy, 3).tolist(),
            'centroid_source_mm': np.round(world, 3).tolist(),
            'inradius_source_mm': None if world_in is None else np.round(world_in, 3).tolist(),
            'polys_xy': polys}


def section_profile(solid, origin, insertion, stations=13, step=INRADIUS_STEP_MM):
    """Cross-sections along the chord plus the summary the design needs."""
    chord = np.asarray(insertion, float) - np.asarray(origin, float)
    length = float(np.linalg.norm(chord))
    rows = []
    for k in range(stations):
        t = 0.02 + 0.96 * k / max(1, stations - 1)
        sec = section_at(solid, origin, insertion, t, step=step)
        if sec is None:
            rows.append({'t': round(t, 4), 'z_mm': round(t * length, 2), 'area_mm2': 0.0,
                         'centroid_xy_mm': None, 'chord_offset_mm': None,
                         'wall_mm': None, 'inradius_mm': None})
            continue
        rows.append({k2: v for k2, v in sec.items() if k2 != 'polys_xy'})
    valid = [r for r in rows if r['wall_mm'] is not None]
    if not valid:
        return {'length_mm': round(length, 2), 'stations': rows, 'sections_inside': 0}
    walls = [r['wall_mm'] for r in valid]
    offs = [r['chord_offset_mm'] for r in valid]
    areas = [r['area_mm2'] for r in valid]
    radii = [r['inradius_mm'] for r in valid if r['inradius_mm'] is not None]
    belly = [r for r in valid if 0.15 <= r['t'] <= 0.85]
    belly_walls = [r['wall_mm'] for r in belly]
    return {'length_mm': round(length, 2),
            'sections_total': len(rows),
            'sections_inside': len(valid),
            'wall_min_mm': round(min(walls), 2),
            'wall_median_mm': round(float(np.median(walls)), 2),
            'wall_min_at_t': min(valid, key=lambda r: r['wall_mm'])['t'],
            'wall_min_belly_mm': round(min(belly_walls), 2) if belly_walls else None,
            'wall_min_belly_at_t': (min(belly, key=lambda r: r['wall_mm'])['t']
                                    if belly else None),
            'inradius_max_mm': round(max(radii), 2) if radii else None,
            'inradius_median_mm': round(float(np.median(radii)), 2) if radii else None,
            'chord_offset_max_mm': round(max(offs), 2),
            'area_max_mm2': round(max(areas), 1),
            'area_median_mm2': round(float(np.median(areas)), 1),
            'stations': rows}
