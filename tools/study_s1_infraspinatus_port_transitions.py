"""Screen smooth port transitions with the existing FJ1500 endpoints fixed."""

from itertools import product

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_smooth_tendon import centreline
from study_s1_infraspinatus_wire_travel import route_metrics
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    wire_file = folder / "FJ1500_wire_travel_and_friction.json"
    wire = read(wire_file)
    deformation = load_reference(wire["deformation"])
    rest = load_reference(deformation["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    prior = load_reference(bias["reference_equilibrium"])
    assembly = load_reference(prior["shoulder_assembly"])
    path = np.array(rest["rest_path_source_mm"])
    lattice = load_component(rest["candidate"])
    keys = [k for k in assembly["sources"] if k.startswith(("bone:", "GH_cap_", "origin:FJ1500/"))
            or k in ("cap:FJ1500/origin", "free_terminal:FJ1500", "stop:FJ1500/distal_stop_blank")]
    solids = {k: load_component(assembly["sources"][k]) for k in keys}
    for key in ("free_terminal:FJ1500", "stop:FJ1500/distal_stop_blank"):
        solids[key] = solids[key].translate(rest["rest_offset_source_mm"])
    original_wire = A.D.outer_envelope(path, .30)
    original_hits = overlaps(original_wire, bounded(solids))
    diagnostics = []
    distal_recut = []
    for label, port, slices in (("proximal", 2, (10, 14, 18, 22, 26)),
                                 ("distal", len(path)-2, (len(path)-11, len(path)-15, len(path)-19, len(path)-23, len(path)-27))):
        # Work in the outward-to-belly direction at both ends.
        oriented = path if label == "proximal" else path[::-1]
        port = port if label == "proximal" else len(path)-1-port
        stations = slices if label == "proximal" else tuple(len(path)-1-i for i in slices)
        targets = bounded(solids)
        candidates = []
        geometry_candidates = []
        for join, (ratio_a, ratio_b), roll in product(stations, product((.25,.5,.8), repeat=2), (0., -12., 12.)):
            start, end = oriented[[port,join]]
            chord = np.linalg.norm(end-start)
            lead = A.D.unit(start-oriented[port-1])
            forward = A.D.unit(oriented[join+1]-end)
            # A rolled interior handle explores clearance while the port
            # station and its actual rigid bore direction are held fixed.
            normal = np.cross(lead, forward)
            if np.linalg.norm(normal)>1e-8:
                forward = Rotation.from_rotvec(A.D.unit(normal)*np.deg2rad(roll)).apply(forward)
            points = centreline([start,start+ratio_a*chord*lead,end-ratio_b*chord*forward,end], 65)
            full = np.vstack([oriented[:port],points,oriented[join+1:]])
            _, turning, radii = route_metrics(full)
            # Include the two seams, whose indices are in the full route.
            transition = slice(port-1, port+64)
            min_radius = float(radii[transition].min())
            row = dict(join_original_index=join if label == "proximal" else len(path)-1-join,
                       handle_ratios=[ratio_a,ratio_b], exit_handle_roll_deg=roll, minimum_transition_radius_mm=min_radius,
                       transition_turn_deg=float(np.rad2deg(turning[transition].sum())))
            candidates.append(row)
            geometry_candidates.append((row,points))
        for row, points in sorted(geometry_candidates, key=lambda p: -p[0]["minimum_transition_radius_mm"])[:15]:
            wire_solid = A.D.outer_envelope(points,.30)
            row["wire_rigid_hits_mm3"] = overlaps(wire_solid,targets)
            row["lattice_stock_requiring_channel_recut_mm3"] = float((wire_solid^lattice).volume())
            row["centreline_source_mm"] = points.tolist()
        diagnostics.append(dict(end=label, held_port_source_mm=oriented[port].tolist(), candidates=candidates))
        if label == "distal":
            mf = A.I._manifold()
            for row, points in sorted(geometry_candidates, key=lambda p: -p[0]["minimum_transition_radius_mm"])[:15]:
                full_path = np.vstack([oriented[:port], points, oriented[len(path)-1-row["join_original_index"]+1:]])
                for wall in (.55, .75):
                    lumen = A.D.outer_envelope(points, .86)
                    new_stock = A.D.outer_envelope(points, .86+wall)
                    # Added material must stay out of nearby rigid parts. Existing
                    # material is tracked separately for the protected stop seat.
                    rigid_neighbourhood = mf.Manifold.batch_boolean(list(solids.values()), mf.OpType.Add)
                    added = new_stock-lattice-rigid_neighbourhood
                    part = ((lattice+added)-lumen).as_original().simplify(1e-6)
                    old_stop = solids["stop:FJ1500/distal_stop_blank"]
                    guard = (old_stop.minkowski_sum(mf.Manifold.cube([1.]*3,True))^lattice)
                    lost = float((guard-part).volume())
                    remaining = float((A.D.outer_envelope(full_path,.30)^part).volume())
                    wall_loss = float(((new_stock-lumen)-part).volume())
                    components = sorted((float(p.volume()) for p in part.decompose() if p.volume()>1e-5), reverse=True)
                    distal_recut.append(dict(join_index=row["join_original_index"], handle_ratios=row["handle_ratios"],
                                             exit_handle_roll_deg=row["exit_handle_roll_deg"], radius_mm=row["minimum_transition_radius_mm"],
                                             nominal_wall_mm=wall, connected_count=len(components), component_volumes_mm3=components,
                                             protected_stop_loss_mm3=lost, nominal_wire_overlap_mm3=remaining,
                                             clipped_wall_volume_mm3=wall_loss, removed_from_TPU_mm3=float((lattice-part).volume())))
    save(folder/"FJ1500_port_transition_screen.json", dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_port_transitions.py"), wire_baseline=provenance(wire_file),
        original_nominal_rigid_hits_mm3=original_hits, ends=diagnostics, distal_lumen_and_wall_trials=distal_recut,
        scope="Static port centre and actual exit tangent held fixed. Cubic transition curvature sampled, including join seams."
              " Best 15 per end checked against original rigid ports, supports, bone and current TPU."
              " Distal recut and added stock are screened as virtual Boolean solids."
              " Stop-seat preservation, wall loss and connectivity are reported; no mesh is adopted."))
    for q in diagnostics:
        checked = [r for r in q["candidates"] if "wire_rigid_hits_mm3" in r]
        print(q["end"],len(q["candidates"]),max(r["minimum_transition_radius_mm"] for r in checked),
              sum(not r["wire_rigid_hits_mm3"] for r in checked),flush=True)


if __name__ == "__main__":
    main()
