"""Compare the old wire with frictionless local wrapping around the sixth rib."""

import heapq

import numpy as np
from scipy.spatial import ConvexHull

from refine_scapular_routing import RIB, support_vertices
from scapular_cable_control import (
    FOLDER, OUTPUT, ROOT, RADIUS_MM, Geometry, length, mesh_geom, posed_points, read, rod, save,
)
from study_elbow_series_package import bounded, overlaps


def cross_section(triangles, origin, normal):
    edges = triangles[:, [1, 2, 0], :] - triangles
    left = (triangles - origin) @ normal
    right = left[:, [1, 2, 0]]
    across = left * right <= 0.
    across &= np.abs(left - right) > 1e-9
    return triangles[across] + edges[across] * (left[across] / (left - right)[across])[:, None]


def enters_polytope(start, end, sides):
    lower, upper = 0., 1.
    for row in sides:
        slope = row[:2] @ (end - start)
        value = row[:2] @ start + row[2] + 1e-7
        if abs(slope) < 1e-12:
            if value >= 0.:
                return False
        elif slope > 0.:
            upper = min(upper, -value / slope)
        else:
            lower = max(lower, -value / slope)
    return upper - lower > 1e-9


def exterior_path(start, end, vertices, sides):
    nodes = np.vstack([start, end, vertices[vertices[:, 1] >= 0.]])
    edges = [[] for _ in nodes]
    for i, a in enumerate(nodes):
        for j in range(i + 1, len(nodes)):
            if not enters_polytope(a, nodes[j], sides):
                span = float(np.linalg.norm(a - nodes[j]))
                edges[i].append((j, span))
                edges[j].append((i, span))
    queue = [(0., 0, (0,))]
    seen = set()
    while queue:
        distance, index, order = heapq.heappop(queue)
        if index == 1:
            return nodes[list(order)]
        if index in seen:
            continue
        seen.add(index)
        for successor, cost in edges[index]:
            if successor not in seen:
                heapq.heappush(queue, (distance + cost, successor, (*order, successor)))
    return None


def wrapped_span(start, end, triangles):
    span = length(np.array([start, end]))
    tangent = (end - start) / span
    lateral = np.array([-1., 0., 0.])
    lateral -= tangent * (lateral @ tangent)
    lateral /= np.linalg.norm(lateral)
    plane = np.stack([tangent, lateral])
    outline = cross_section(triangles, start, np.cross(tangent, lateral))
    if len(outline) < 3:
        return np.array([start, end])
    projected = np.unique(np.round((outline - start) @ plane.T, 7), axis=0)
    if len(projected) < 3 or np.linalg.matrix_rank(projected - projected.mean(axis=0)) < 2:
        return np.array([start, end])
    hull = ConvexHull(projected)
    polygon = projected[hull.vertices]
    edges = np.roll(polygon, -1, axis=0) - polygon
    normals = np.column_stack([edges[:, 1], -edges[:, 0]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    offsets = np.einsum("ij,ij->i", polygon, normals) + RADIUS_MM + .02
    expanded = np.array([np.linalg.solve(np.stack([normals[i - 1], normals[i]]), [offsets[i - 1], offsets[i]])
                         for i in range(len(polygon))])
    path = exterior_path(np.zeros(2), np.array([span, 0.]), expanded, np.column_stack([normals, -offsets]))
    return None if path is None else start + path @ plane


def main():
    geo = Geometry()
    source = read(OUTPUT / "scapular_cable_control.json")
    serratus = source["routes"]["serratus_middle"]
    action = read(FOLDER / "whole_arm_anatomical_actions.json")["scapular"]["upward_rotation"]
    entry = geo.sources[RIB]
    _, faces = mesh_geom.load_obj(ROOT / entry["file"])
    triangles = support_vertices(geo, RIB)[faces]
    # This patch represents the local obstacle encountered by the old final span.
    box_center = np.array([-98., -31., 1231.])
    triangles = triangles[np.all(abs(triangles.mean(axis=1) - box_center) <= [29., 25., 23.], axis=1)]
    rows = []
    for angle in np.arange(0., 10.1, .5):
        data = geo.prescribed(action, angle)
        first, start, end = posed_points(data, geo.origins, serratus)
        path = wrapped_span(start, end, triangles)
        if path is None:
            rows.append(dict(angle_deg=angle, status="no lateral path in section"))
            continue
        rigid, wires = geo.posed_solids(data, excluded={"FJ1459_g3"})
        obstacles = bounded({**rigid, **wires})
        geo.cached.clear()
        hits = {}
        for a, b in zip(path, path[1:]):
            for obstacle, volume in overlaps(rod(a, b, RADIUS_MM), obstacles).items():
                hits[obstacle] = hits.get(obstacle, 0.) + volume
        rows.append(dict(angle_deg=angle, status="screened", point_count=len(path),
                         final_span_points_mm=path.tolist(), total_length_mm=length(np.vstack([first, path])),
                         hits=hits))
    save(OUTPUT / "ideal_rib_wrap_probe.json",
         dict(method="Frictionless taut path around the expanded convex section of the local sixth-rib patch; "
                     "exterior branch only, subject to a separate 3D rigid clearance check.",
              validated_motion=False, rows=rows))
    for row in rows:
        print(row["angle_deg"], row["status"], row.get("point_count"), list(row.get("hits", {})))


if __name__ == "__main__":
    main()
