"""Compare tension-only opposed TPU bands as a physical source of coupling."""

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_y_retained_lead_clearance import rigid_transform
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from check_shoulder_front_liner_keepers import is_compliant
from study_s1_infraspinatus_cap_relief import HEAD
from study_s1_infraspinatus_tpu_bias import collision_state, rest_geometry
from study_s1_local_force_geometry import load_reference
from study_s1_y_terminal_equilibrium import solve_branch


def band_model(axial, sideways, span, preload):
    axis = (sideways - axial) / np.sqrt(2.)
    points = np.array([-axis, axis]) * span
    stiffness = .2 / (1. - preload)
    in_plane = stiffness * (1. + preload)
    base = np.array([.3 - in_plane, .4 - in_plane, .4 - 2. * stiffness * preload])
    if np.any(base < 0.):
        raise ValueError("Matching this preload would require negative base stiffness")
    frame = np.stack([axial, sideways, np.cross(axial, sideways)])
    rest_length = span * (1. - preload)

    def loads(displacement):
        return stiffness * np.maximum(np.linalg.norm(points - displacement, axis=1) - rest_length, 0.)

    def energy(displacement):
        axial_coordinates = frame @ displacement
        return .5 * (base @ (axial_coordinates**2) + (loads(displacement) @ loads(displacement)) / stiffness)

    return points, stiffness, base, rest_length, loads, energy


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    coupling_path = folder / "FJ1500_tpu_coupling_screen.json"
    coupling = read(coupling_path)
    path = load_reference(coupling["prior_path"])
    rest = load_reference(path["rest_geometry"])
    bias = load_reference(rest["bias_screen"])
    previous = load_reference(bias["reference_equilibrium"])
    geometry = load_reference(previous["local_geometry"])
    assembly = load_reference(previous["shoulder_assembly"])
    model = next(row for row in load_reference(previous["kinematics"])["models"] if row["head"] == HEAD)
    ports = next(row["ports"] for row in load_reference(geometry["port_geometry"])["routes"] if row["head"] == HEAD)
    sleeve = load_reference(path["sleeve"])
    sources = dict(assembly["sources"])
    key = f"link:{HEAD}/muscle_uncrimped_sleeve"
    if sources[key] != sleeve["original"]:
        raise ValueError("Compact sleeve and S1 assembly differ")
    sources[key] = sleeve["selected"]["source"]
    solids = {key: load_component(source) for key, source in sources.items()
              if not is_compliant(key) and not key.startswith(("bench:", "fixture:", "clamp:", "clavicle:"))
              and not key.endswith("/tendon_wire")}
    centre = np.array(geometry["gh_centre_source_mm"])
    flexion = np.array(geometry["axes_source"][geometry["axes"].index("flexion")])
    pin_a, pin_b = (np.array(model[f"neutral_{end}_pin_source_mm"]) for end in ("cap", "bone"))
    axis_a, axis_b = (np.array(model[f"{end}_pin_axis_source"]) for end in ("cap", "bone"))
    offset = np.array(rest["rest_offset_source_mm"])
    adjusted, adjusted_ports = rest_geometry(model, ports, offset)
    origin, stop = [np.array(adjusted_ports[port][face]) for port, face in
                    (("origin", "inner_face_source_mm"), ("distal", "outer_face_source_mm"))]
    axial = (origin - stop) / np.linalg.norm(origin - stop)
    sideways = offset - (offset @ axial) * axial
    sideways /= np.linalg.norm(sideways)
    owners = assembly["intended_body_owners"]
    cases = []
    angles = (0., 30., 60., 64., 68., 75., 85., 96.)
    for span, preload in ((12., .08), (18., .08), (24., .08), (18., .16), (24., .16)):
        points, stiffness, base, rest_length, loads, energy = band_model(axial, sideways, span, preload)
        anchors = points + pin_a + offset
        target = np.array([[.3, -.2, 0.], [-.2, .4, 0.], [0., 0., .4]])
        frame = np.stack([axial, sideways, np.cross(axial, sideways)])
        probe = np.zeros((3, 3))
        for i, a in enumerate(frame * .001):
            for j, b in enumerate(frame * .001):
                probe[i, j] = (energy(a + b) - energy(a - b) - energy(-a + b) + energy(-a - b)) / (4e-6)
        if not np.allclose(probe, target, atol=1e-7, rtol=0.):
            raise ValueError("Band tangent response did not match the coupling target")
        records = []
        for angle in angles:
            gh = Rotation.from_rotvec(flexion * np.deg2rad(angle))
            bone_frame = rigid_transform(gh, centre, centre)
            states = solve_branch(adjusted, adjusted_ports, centre, gh, bias["stiffness"], [0., 5., 10.],
                                  retain_leads=True, translation_energy=energy)
            for state in states:
                pin, hits = collision_state(solids, owners, state, offset, gh, centre, pin_a, pin_b, axis_a, axis_b)
                force = loads(np.array(state["cap_translation_source_mm"]))
                band_hits = []
                for band, (anchor, tension) in enumerate(zip(anchors, force)):
                    if tension <= 1e-5:
                        continue
                    # 0.7 mm is a nominal clearance probe for a tension band,
                    # not a selected material section or a modeled sleeve.
                    sweep = rod(anchor, pin, .7)
                    for name in ("bone:scapula", "bone:humerus", "GH_cap_0", "GH_cap_1"):
                        part = solids[name]
                        if owners[name] == "humerus_ideal_attachment":
                            part = part.transform(bone_frame)
                        overlap = float(A.I.intersection_volume(sweep, part))
                        if overlap > .001:
                            band_hits.append(dict(band=band, part=name, overlap_mm3=overlap))
                records.append(dict(GH_deg=angle, input_N=state["commanded_muscle_tension_N"],
                                    pin_source_mm=pin.tolist(), cap_and_eye_hits=hits, band_hits=band_hits,
                                    band_forces_N=force.tolist(),
                                    equilibrium_residual_N=state["stationarity_residual_N"]))
        result = dict(span_mm=span, preload_ratio=preload, band_stiffness_N_per_mm=stiffness,
                      complementary_TPU_stiffness_N_per_mm=base.tolist(),
                      neutral_preload_per_band_N=stiffness * (span - rest_length),
                      rest_length_mm=rest_length, fixed_anchor_candidates_source_mm=anchors.tolist(),
                      tangent_stiffness_N_per_mm=probe.tolist(),
                      conflicting_rigid_states=sum(bool(row["cap_and_eye_hits"]) for row in records),
                      conflicting_band_states=sum(bool(row["band_hits"]) for row in records), rows=records)
        cases.append(result)
        print(span, preload, result["conflicting_rigid_states"], result["conflicting_band_states"], flush=True)
    save(folder / "FJ1500_tension_band_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_infraspinatus_tension_bands.py"),
        solver=provenance(ROOT / "tools/study_s1_y_terminal_equilibrium.py"),
        coupling=provenance(coupling_path), cases=cases,
        scope="Opposed tension-only bands with equilibrated neutral preload, plus complementary TPU stiffness."
              " Fixed support candidates are relative to scapula but no bracket is designed."
              " Straight tension bands are screened as 0.7 mm-radius cylinders only when taut."
              " Terminal strap attachment, fatigue, wire curvature and integrated dynamics are pending."))


if __name__ == "__main__":
    main()
