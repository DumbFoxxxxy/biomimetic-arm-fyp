"""Check the central strand and measure passive support requirements for FJ1500."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_elbow_series_package import bounded, overlaps
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_smooth_tendon import centreline
from study_s1_infraspinatus_tpu_bias import rest_geometry
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_outward_path_screen.json"
    report = read(path_file)
    smooth = load_reference(report["parent"])
    coupling = load_reference(smooth["coupling"])
    path = load_reference(coupling["prior_path"])
    rest = load_reference(path["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    prior = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(prior["local_geometry"])
    assembly = load_reference(report["shoulder_assembly"])
    model = next(row for row in load_reference(prior["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(report["input_route"])["routes"] if row["head"] == HEAD)
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sleeve = load_reference(path["sleeve"])
    if sleeve["original"] != assembly["sources"][sleeve_key]:
        raise ValueError("Compact sleeve no longer matches the shoulder")
    solids = {key: load_component(source) for key, source in {**assembly["sources"], sleeve_key: sleeve["selected"]["source"]}.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    owners = assembly["intended_body_owners"]
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pins = [np.array(model[f"neutral_{k}_pin_source_mm"]) for k in ("cap", "bone")]
    axes = [np.array(model[f"{k}_pin_axis_source"]) for k in ("cap", "bone")]
    outward = np.array(report["direction_source"])
    chosen = report["selected"]
    _, ports = rest_geometry(model, route["ports"], outward*chosen["rest_mm"])
    origin, stop = (np.array(ports[key][face]) for key, face in (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    axial = A.D.unit(origin-stop)
    raw_origin, raw_stop = (np.array(route["ports"][key][face])
                            for key, face in (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm")))
    raw_axis = A.D.unit(raw_origin-raw_stop)
    direction = Rotation.from_rotvec(raw_axis*np.deg2rad(chosen["coupling_roll_deg"])).apply(outward)
    lateral = A.D.unit(direction - (direction @ axial)*axial)
    other = np.cross(axial, lateral)
    stiffness_axial, stiffness_lateral = bias["stiffness"][:2]
    gain = chosen["gain"]
    stiffness = np.array([[stiffness_axial + gain**2*stiffness_lateral, -gain*stiffness_lateral, 0.],
                          [-gain*stiffness_lateral, stiffness_lateral, 0.], [0., 0., stiffness_lateral]])
    records = []
    for branch, states in (("outward", report["outward"]), ("return", report["return_rows"])):
        for state in states:
            angle = state["GH_deg"]
            load = state["commanded_muscle_tension_N"]
            limit = report["limits"][f"{load:g}_{branch}"]["last_clear_deg"]
            if limit is None or angle > limit:
                continue
            gh = Rotation.from_rotvec(flexion*np.deg2rad(angle))
            cap = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
            pin = np.array(state["cap_pin_source_mm"])
            turn_a, turn_b = np.deg2rad(state["swivel_angles_deg"])
            body = rigid_transform(gh, centre, centre)
            cap_pose = rigid_transform(cap, pins[0], pin)
            eye_a = rigid_transform(cap*Rotation.from_rotvec(axes[0]*turn_a), pins[0], pin)
            eye_b = rigid_transform(gh*Rotation.from_rotvec(axes[1]*turn_b), pins[1], centre+gh.apply(pins[1]-centre))
            exits, tips = (np.array(state[key]) for key in ("tendon_exits_source_mm", "lead_tips_source_mm"))
            scale = np.linalg.norm(tips[1]-tips[0])/3.
            points = centreline([tips[0], tips[0]+scale*A.D.unit(tips[0]-exits[0]),
                                 tips[1]+scale*A.D.unit(tips[1]-exits[1]), tips[1]], 65)
            wire = A.D.outer_envelope(points, .30)
            targets = {}
            for key, part in solids.items():
                if key in (f"link:{HEAD}/muscle_groove_ring", sleeve_key):
                    part = part.transform(eye_a)
                elif key in (f"link:{HEAD}/bone_groove_ring", f"link:{HEAD}/bone_uncrimped_sleeve"):
                    part = part.transform(eye_b)
                elif owners[key] == f"free_terminal:{HEAD}":
                    part = part.transform(cap_pose)
                elif owners[key] == "humerus_ideal_attachment":
                    part = part.transform(body)
                targets[key] = part
            displacement = np.stack([axial, lateral, other]) @ state["cap_translation_source_mm"]
            records.append(dict(branch=branch, GH_deg=angle, input_N=load,
                                local_cap_displacement_mm=displacement.tolist(),
                                translational_TPU_reaction_N=(-stiffness@displacement).tolist(),
                                short_tendon_tension_N=state["short_tendon_tension_N"],
                                pull_chord_change_mm=state["muscle_proxy_chord_change_mm"],
                                central_radius_mm=state["minimum_sampled_radius_mm"],
                                central_wire_hits_mm3=overlaps(wire, bounded(targets))))
        print(branch, len(records), flush=True)
    save(folder / "FJ1500_outward_wire_and_support.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_infraspinatus_outward_wire.py"), trajectory=provenance(path_file),
        support_basis_source=[axial.tolist(), lateral.tolist(), other.tolist()], stiffness_N_per_mm=stiffness.tolist(),
        records=records, scope="Central 0.60 mm strand; two retained crimp lead zones are outside the wire-overlap check."
              " Translational reactions are calculated from the assumed passive law."
              " Wire friction, actual stiffness, cap support geometry and free driving are unverified."))


if __name__ == "__main__":
    main()
