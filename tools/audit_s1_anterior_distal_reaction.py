"""Locate solid stock in the FJ1468 distal turn load direction."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import load_component, save
from channel_first_bone_contact import BoneContactSurface
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    audit_ref = provenance(folder / "FJ1468_anchor_path_audit.json")
    audit = load_reference(audit_ref)
    assembly = load_reference(audit["anchor_register"])["neutral_rigid_sources"]
    route = load_reference(audit["revised_route"])
    assembly["tpu:existing_front_guide"] = route["lined_tpu"]
    points = np.asarray(route["neutral_path_source_mm"])
    nodes = {row["peak_turn_node"] for row in audit["cases"]}
    if len(nodes) != 1 or any(row["peak_turn_distal_blend_weight"] != 1. for row in audit["cases"]):
        raise ValueError("The peak turn does not consistently follow the humerus")
    node = nodes.pop()
    center = points[node]
    tangent = np.diff(points[node - 1:node + 2], axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    load = tangent[1] - tangent[0]
    magnitude = float(np.linalg.norm(load))
    direction = load / magnitude
    side = np.cross(tangent.sum(axis=0), direction)
    side /= np.linalg.norm(side)
    nearby = ("cap:FJ1468/distal", "stop:FJ1468_stop_blank", "bone:humerus",
              "elbow:FJ1486_bone_adapter", "tpu:existing_front_guide")
    rows = {}
    for key in nearby:
        solid = load_component(assembly[key])
        surface = BoneContactSurface(solid, band_mm=15.)
        distance, _ = surface.nearest(center)
        rays = []
        for theta in (-25., -10., 0., 10., 25.):
            angle = np.deg2rad(theta)
            radial = direction * np.cos(angle) + side * np.sin(angle)
            hits = solid.ray_cast(center, center + 12. * radial)
            entering = [hit for hit in hits if np.dot(hit.normal, radial) < 0.]
            station = float((np.asarray(entering[0].position) - center) @ radial) if entering else None
            rays.append(dict(offset_deg=theta, entry_from_wire_axis_mm=station,
                             gap_from_nominal_wire_mm=station - .30 if station is not None else None))
        rows[key] = dict(signed_nearest_distance_mm=float(distance), rays=rays, source=assembly[key])
    save(folder / "FJ1468_distal_reaction_support.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_anterior_distal_reaction.py"), route_audit=audit_ref,
        node=node, point_source_mm=center.tolist(), reaction_per_tension=magnitude,
        force_on_support_direction_source=direction.tolist(), components=rows,
        scope="Neutral geometric rays in the equal-tension turn load direction and a transverse fan. "
              "No contact equilibrium, material stiffness, or fastener capacity inferred. "
              "Nearest distances are capped at 15 mm. The existing TPU mesh is checked as evidence "
              "of the present geometry; it is not changed.",
    ))
    print({key: (row["signed_nearest_distance_mm"], row["rays"][2]["entry_from_wire_axis_mm"])
           for key, row in rows.items()}, flush=True)


if __name__ == "__main__":
    main()
