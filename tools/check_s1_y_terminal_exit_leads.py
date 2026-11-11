"""Check whether solved chord states leave room for the actual crimp exit leads."""

import json

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load"
    equilibrium_file = folder / "y_terminal_equilibrium.json"
    equilibrium = read(equilibrium_file)
    geometry = load_reference(equilibrium["local_geometry"])
    models = {row["head"]: row for row in load_reference(equilibrium["kinematics"])["models"]}
    axes = dict(zip(geometry["axes"], np.array(geometry["axes_source"])))
    rows = []
    for case in equilibrium["cases"]:
        model = models[case["head"]]
        path = np.array(model["active_strand_centreline_source_mm"])
        lead_a, lead_b = path[1] - path[0], path[-2] - path[-1]
        lead_total = np.linalg.norm(lead_a) + np.linalg.norm(lead_b)
        axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
        gh = Rotation.from_rotvec(axes.get(case["prescribed_axis"], axes["flexion"])
                                  * np.deg2rad(case["prescribed_GH_deg"]))
        for step, state in enumerate(case["rows"]):
            point_a, point_b = np.array(state["tendon_exits_source_mm"])
            cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
            eye_a, eye_b = np.deg2rad(state["swivel_angles_deg"])
            tip_a = point_a + cap.apply(Rotation.from_rotvec(axis_a * eye_a).apply(lead_a))
            tip_b = point_b + gh.apply(Rotation.from_rotvec(axis_b * eye_b).apply(lead_b))
            lower_bound = lead_total + np.linalg.norm(tip_a - tip_b)
            missing = lower_bound - model["nominal_free_strand_material_length_mm"]
            rows.append(dict(head=case["head"], prescribed_axis=case["prescribed_axis"],
                             prescribed_GH_deg=case["prescribed_GH_deg"], stiffness=case["stiffness"],
                             step=step, commanded_muscle_tension_N=state["commanded_muscle_tension_N"],
                             short_tendon_tension_N=state["short_tendon_tension_N"],
                             lead_constrained_lower_bound_mm=lower_bound,
                             additional_length_required_at_least_mm=float(max(0., missing)),
                             lower_bound_exceeds_assumed_installed_length=bool(missing > 1e-5)))
    loaded = [row for row in rows if row["short_tendon_tension_N"] > 1e-5]
    conflicts = [row for row in loaded if row["lower_bound_exceeds_assumed_installed_length"]]
    report = dict(generator=provenance(ROOT / "tools/check_s1_y_terminal_exit_leads.py"),
                  equilibrium=provenance(equilibrium_file), rows=rows,
                  loaded_states=len(loaded), loaded_length_conflicts=len(conflicts),
                  scope="Each source straight exit lead is retained. The distance between lead tips plus lead"
                        " lengths is a necessary lower bound. Success is not tangent continuity or clearance."
                        " Excess length applies to the prescribed equilibrium state and assumed installed length.")
    (folder / "y_terminal_exit_lead_check.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{len(conflicts)} of {len(loaded)} tensioned states fail the retained-lead length lower bound")


if __name__ == "__main__":
    main()
