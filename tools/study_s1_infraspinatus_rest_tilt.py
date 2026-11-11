"""Screen a common cap rest tilt for the laterally offset FJ1500 TPU."""

import numpy as np
from scipy.spatial.transform import Rotation

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def tilt_geometry(model, ports, offset, tilt):
    adjusted, adjusted_ports = rest_geometry(model, ports, offset)
    pin = np.array(adjusted["neutral_cap_pin_source_mm"])
    key = "neutral_crimp_outer_exit_A_source_mm"
    adjusted[key] = (pin + tilt.apply(np.array(adjusted[key]) - pin)).tolist()
    adjusted["cap_pin_axis_source"] = tilt.apply(adjusted["cap_pin_axis_source"]).tolist()
    path = np.array(adjusted["active_strand_centreline_source_mm"])
    path[:2] = pin + tilt.apply(path[:2] - pin)
    adjusted["active_strand_centreline_source_mm"] = path.tolist()
    key = "outer_face_source_mm"
    adjusted_ports["distal"][key] = (pin + tilt.apply(np.array(adjusted_ports["distal"][key]) - pin)).tolist()
    return adjusted, adjusted_ports


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    path_file = folder / "FJ1500_flexion_path.json"
    path = read(path_file)
    rest = load_reference(path["rest_geometry"])
    screen = load_reference(rest["bias_screen"])
    previous = load_reference(screen["reference_equilibrium"])
    geometry = load_reference(previous["local_geometry"])
    assembly = load_reference(previous["shoulder_assembly"])
    model = next(row for row in load_reference(previous["kinematics"])["models"] if row["head"] == HEAD)
    ports = next(row["ports"] for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sources = dict(assembly["sources"])
    sleeve_key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    sleeve = load_reference(path["sleeve"])
    if sources[sleeve_key] != sleeve["original"]:
        raise ValueError("Compact sleeve belongs to a different assembly")
    sources[sleeve_key] = sleeve["selected"]["source"]
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    offset = np.array(rest["rest_offset_source_mm"])
    stop = np.array(ports["distal"]["outer_face_source_mm"])
    tilt_axis = np.cross(stop - pin_a, offset)
    tilt_axis /= np.linalg.norm(tilt_axis)
    rows = []
    for angle in (0., -15., 10., 20., 30.):
        tilt = Rotation.from_rotvec(tilt_axis * np.deg2rad(angle))
        adjusted, route = tilt_geometry(model, ports, offset, tilt)
        for gh_deg in (0., 60., 64., 68., 70., 75., 80., 96.):
            gh = Rotation.from_rotvec(flexion * np.deg2rad(gh_deg))
            for state in solve_branch(adjusted, route, centre, gh, screen["stiffness"], [0., 5., 10.], retain_leads=True):
                delta = Rotation.from_rotvec(np.deg2rad(state["cap_rotation_vector_deg"]))
                state["cap_rotation_vector_deg"] = np.rad2deg((delta * tilt).as_rotvec()).tolist()
                pin, hits = collision_state(solids, assembly["intended_body_owners"], state, offset, gh,
                                            centre, pin_a, pin_b, axis_a, axis_b)
                rows.append(dict(rest_tilt_deg=angle, GH_flexion_deg=gh_deg, cap_pin_source_mm=pin.tolist(),
                                 hits=hits, **state))
        tested = [row for row in rows if row["rest_tilt_deg"] == angle and row["commanded_muscle_tension_N"] == 10.]
        print("Tilt", angle, "GH and rigid hits", [(r["GH_flexion_deg"], len(r["hits"])) for r in tested], flush=True)
    save(folder / "FJ1500_rest_tilt_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_rest_tilt.py"), reference=provenance(path_file),
        common_rest_tilt_axis_source=tilt_axis.tolist(), common_lateral_offset_source_mm=offset.tolist(), rows=rows,
        scope="Prescribed GH and solved terminal equilibrium at several probes. The rest tilt rotates the cap,"
              " actual stop and muscle-side eye geometry consistently. No TPU lattice export for the tilt."
              " The prior short-tendon straight-tip lower bound and diagnostic rigid-only checks remain."))


if __name__ == "__main__":
    main()
