"""Compare shared passive cap trajectories with the smooth FJ1500 short strand."""

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    parent_file = folder / "FJ1500_smooth_tendon_screen.json"
    parent = read(parent_file)
    coupling = load_reference(parent["coupling"])
    path = load_reference(coupling["prior_path"])
    rest = load_reference(path["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    prior = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(prior["local_geometry"])
    assembly = load_reference(prior["shoulder_assembly"])
    model = next(row for row in load_reference(prior["kinematics"])["models"] if row["head"] == HEAD)
    route = next(row for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sleeve = load_reference(path["sleeve"])
    key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    if assembly["sources"][key] != sleeve["original"]:
        raise ValueError("Sleeve and assembly mismatch")
    sources = {**assembly["sources"], key: sleeve["selected"]["source"]}
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    outward = np.array(coupling["coupling_direction_source"])
    points = [np.array(route["ports"][k][face]) for k, face in (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm"))]
    axial = points[0]-points[1]
    axial /= np.linalg.norm(axial)
    pins = [np.array(model[f"neutral_{k}_pin_source_mm"]) for k in ("cap", "bone")]
    pin_axes = [np.array(model[f"{k}_pin_axis_source"]) for k in ("cap", "bone")]

    def evaluate(rest_mm, roll_deg, gain, angles, levels):
        offset = outward * rest_mm
        direction = Rotation.from_rotvec(axial*np.deg2rad(roll_deg)).apply(outward)
        adjusted, ports = rest_geometry(model, route["ports"], offset)
        rows = []
        for angle in angles:
            gh = Rotation.from_rotvec(flexion*np.deg2rad(angle))
            states = solve_branch(adjusted, ports, centre, gh, bias["stiffness"], levels, retain_leads=True,
                                  lateral_coupling=(direction, gain), strand_bending_Nmm2=10.)
            for state in states:
                pin, hits = collision_state(solids, assembly["intended_body_owners"], state, offset, gh, centre, *pins, *pin_axes)
                rows.append(dict(GH_deg=angle, cap_pin_source_mm=pin.tolist(), hits=hits, **state))
        return rows

    trials = []
    cases = [(6., roll, gain) for roll in (0., -30., 30., -60., 60., -90., 90.) for gain in (.5, 1.)]
    cases += [(distance, roll, gain) for distance in (8., 10.) for roll in (0., -45., 45.) for gain in (.5, 1.)]
    for rest_mm, roll, gain in cases:
        rows = evaluate(rest_mm, roll, gain, (0., 30., 60., 68., 75., 85., 96.), [0., 5., 10.])
        hits = [row for row in rows if row["hits"]]
        if (rest_mm, roll, gain) == (6., 0., .5):
            reference = {(r["GH_deg"], r["commanded_muscle_tension_N"]): r for r in parent["arc_equilibrium"]
                         if r["assumed_bending_Nmm2"] == 10.}
            for row in rows:
                peer = reference.get((row["GH_deg"], row["commanded_muscle_tension_N"]))
                if peer and (not np.allclose(row["cap_pin_source_mm"], peer["cap_pin_source_mm"], atol=2e-6, rtol=0.)
                             or bool(row["hits"]) != bool(peer["rigid_terminal_hits"])):
                    raise ValueError("Baseline cap trajectory changed")
        trials.append(dict(rest_mm=rest_mm, coupling_roll_deg=roll, gain=gain, failed_states=len(hits),
                           first_failed_probe_deg=hits[0]["GH_deg"] if hits else None, rows=rows))
        print(rest_mm, roll, gain, 'conflicts', len(hits), flush=True)
    selected = min(trials, key=lambda row: (-(row["first_failed_probe_deg"] if row["first_failed_probe_deg"] is not None else 180.),
                                             row["failed_states"], row["rest_mm"], row["gain"], abs(row["coupling_roll_deg"])))
    args = selected["rest_mm"], selected["coupling_roll_deg"], selected["gain"]
    onward = evaluate(*args, list(range(97)), [0., 5., 10.])
    returning = evaluate(*args, list(range(96, -1, -1)), [10., 5., 0.])
    limits = {}
    for load in (0., 5., 10.):
        for branch, rows in (("outward", onward), ("return", returning)):
            matches = {q["GH_deg"]: q for q in rows if q["commanded_muscle_tension_N"] == load}
            failed = next((angle for angle in range(97) if matches[angle]["hits"]), None)
            detail = []
            if failed is not None and failed > 0:
                angles = [round(failed-1+i*.1, 1) for i in range(1, 10)]
                if branch == "return":
                    angles.reverse()
                detail = evaluate(*args, angles, [0., 5., 10.] if branch == "outward" else [10., 5., 0.])
                matches.update({q["GH_deg"]: q for q in detail if q["commanded_muscle_tension_N"] == load})
            ordered = sorted(matches)
            idx = next((i for i, angle in enumerate(ordered) if matches[angle]["hits"]), len(ordered))
            limits[f"{load:g}_{branch}"] = dict(last_clear_deg=ordered[idx-1] if idx else None,
                                                  first_contact=matches[ordered[idx]] if idx<len(ordered) else None,
                                                  fine=detail)
    save(folder / "FJ1500_outward_path_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_outward_path.py"), parent=provenance(parent_file),
        shoulder_assembly=prior["shoulder_assembly"], input_route=geometry["port_geometry"],
        direction_source=outward.tolist(), trials=trials, selected={k: v for k, v in selected.items() if k != "rows"},
        limits=limits, outward=onward, return_rows=returning,
        scope="Rest offset and passive lateral coupling are shared across poses. Prescribed GH angle."
              " Rigid cap and eye collisions only. The candidate stiffness law is not linked to a built flexure."
              " TPU deformation, actual wire contact, whole arm loading and manufacture remain open."))
    print('selected',args, 'limits', {k:v['last_clear_deg'] for k,v in limits.items()}, flush=True)


if __name__ == "__main__":
    main()
