"""Check the retained return spring and modified hardware on the shoulder–elbow trace."""

import argparse
from itertools import combinations

import mujoco
import numpy as np

from audit_elbow_slider_short_guides import displacement_map, weigh
from audit_s1_series_elbow_integration import MOBILE, find_hits, transforms
from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_series_elbow_return_spring import spring_path, sweep
from build_s1_teres_origin_mount import A, load_component, save
from refine_elbow_short_guide_relief import indexed_solid
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def load_candidate(folder):
    retention_file = folder / "return_spring_pin_retention.json"
    retention = read(retention_file)
    if (not retention["cross_body_slider_sample_count"] or retention["cross_body_slider_hits"]
            or retention["spring_full_stroke_hits"] or any(retention["combined_preassembly_hit_offsets_mm"].values())):
        raise ValueError("Missing locally checked spring retention assembly")
    seats = load_reference(retention["parent"])
    spring = load_reference(seats["parent"])
    size = load_reference(spring["size"])
    prior = load_reference(size["geometry"])
    extended = load_reference(prior["extended_drive"])
    assembly = load_reference(load_reference(extended["budget"])["assembly"])
    mount = load_reference(size["spring_mount"])["spring_mount"]
    sources = {**assembly["sources"], **extended["shafts"]}
    owners = dict(assembly["owners"])
    revisions = {}
    for label, key in (("fixed", "elbow:FJ1477_bone_adapter"), ("moving", "elbow:FJ1477_proximal_rigid_cap")):
        meshes = retention["results"][label]["selected"]["meshes"]
        revisions[key] = meshes["lug"]
        revisions[f"series:{label}_pin"] = meshes["pin"]
        screw = f"series:{label}_dog_screw"
        revisions[screw] = meshes["dog_screw"]
        owners[screw] = owners[key]
    return prior, extended, mount, spring, sources, owners, revisions, provenance(retention_file)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    spool_cases = ("mu_30_original", "mu_30_adjusted", "seating_k5", "seating_k20", "seating_k50")
    mode.add_argument("--anterior-friction", action="store_true",
                        help="Replay the extended-rail 2 N command / mu=0.3 trajectory")
    mode.add_argument("--nominal-rate-rigid", action="store_true",
                      help="Check the nominal sized spring, zero anterior command, against the rigid assembly")
    mode.add_argument("--nominal-rate-flexible", action="store_true",
                      help="Replay the coil, wires and TPU on the nominal sized-spring trajectory")
    mode.add_argument("--spool-rigid", choices=spool_cases,
                      help="Audit a 0.3 friction spool trajectory against the rigid assembly")
    mode.add_argument("--spool-flexible", choices=spool_cases,
                      help="Replay the spring and elbow muscle geometry on a rigid-checked spool trace")
    args = parser.parse_args()
    folder = DEST.parent / "integration_review/s1_series_elbow"
    prior, extended, mount, spring, sources, owners, revisions, retention = load_candidate(folder)
    trace = checked_trace(prior["trace"])
    if (spring["trajectory_frames"] != len(trace) or prior["trace"] != extended["trajectory"]
            or len(prior["records"]) != len(trace)):
        raise ValueError("Inconsistent spring and wire trajectories")
    trace_ref = prior["trace"]
    friction_ref = None
    rate_ref = None
    spool_case = args.spool_rigid or args.spool_flexible
    spool_ref = None
    spool_rigid_ref = None
    if args.anterior_friction:
        friction_ref = provenance(folder / "extended_anterior_friction_integration.json")
        friction = load_reference(friction_ref)
        if friction["model"] != extended["model"]:
            raise ValueError("Friction trajectory uses a different drive model")
        case = next(row for row in friction["cases"] if row["mu"] == .3)
        if case["first_rigid_collision"] or case["first_pin_rib_collision"]:
            raise ValueError("Selected friction trajectory has a known rigid collision")
        trace_ref = case["trace"]
        alternative = checked_trace(trace_ref)
        if [row["time_s"] for row in alternative] != [row["time_s"] for row in trace]:
            raise ValueError("Friction trajectory does not match the baseline sampling grid")
        trace = alternative
    if spool_case:
        study = "seating" if spool_case.startswith("seating_") else "friction"
        spool_ref = provenance(folder / f"nominal_spring_spool_{study}.json")
        spool = load_reference(spool_ref)
        if spool["drive_model"] != extended["model"]:
            raise ValueError("Spool trajectory uses a different drive model")
        case = next(row for row in spool["cases"] if row["case"] == spool_case)
        trace_ref = case["trace"]
        alternative = checked_trace(trace_ref)
        if [row["time_s"] for row in alternative] != [row["time_s"] for row in trace]:
            raise ValueError("Spool trajectory does not match the baseline sampling grid")
        trace = alternative
        if args.spool_flexible:
            spool_rigid_ref = provenance(folder / f"nominal_spring_spool_{spool_case}_rigid_audit.json")
            audit = load_reference(spool_rigid_ref)
            if (audit["trace"] != trace_ref or audit["retention"] != retention or audit["findings"]
                    or audit["spool_study"] != spool_ref
                    or audit["sampled_frames"] != len(trace)):
                raise ValueError("Missing matching clear rigid audit for the spool trajectory")
    if args.nominal_rate_rigid or args.nominal_rate_flexible:
        rate_ref = provenance(folder / "return_spring_rate_dynamics.json")
        rate = load_reference(rate_ref)
        if rate["model"] != extended["model"] or rate["retention"] != retention:
            raise ValueError("Spring-rate trajectory uses different hardware or model")
        nominal_modulus = load_reference(spring["size"])["assumed_shear_modulus_MPa"][1]
        case = next(row for row in rate["cases"] if row["modulus_MPa"] == nominal_modulus
                    and row["anterior_command_peak_N"] == 0.)
        trace_ref = case["trace"]
        alternative = checked_trace(trace_ref)
        if [row["time_s"] for row in alternative] != [row["time_s"] for row in trace]:
            raise ValueError("Spring-rate trajectory does not match the baseline sampling grid")
        trace = alternative
    model_file = ROOT / extended["model"]["file"]
    if provenance(model_file) != extended["model"]:
        raise ValueError("Drive model changed")
    model = mujoco.MjModel.from_xml_string(model_file.read_text(encoding="utf-8"))
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    neutral = {name: data.body(name).xpos.copy() * 1000. for name in MOBILE}
    original = {key: load_component(sources[key]) for key in revisions if key in sources}
    hardware = {key: load_component(source) for key, source in {**sources, **revisions}.items()}
    if args.nominal_rate_rigid or args.spool_rigid:
        keys = list(hardware)
        pairs = np.array([(i, j) for i, j in combinations(range(len(keys)), 2)
                          if owners[keys[i]] != owners[keys[j]]
                          and (owners[keys[i]] in MOBILE or owners[keys[j]] in MOBILE)])
        findings = []
        clear_peak = None
        for index, row in enumerate(trace):
            data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
            mujoco.mj_kinematics(model, data)
            frames = transforms(model, data, neutral)
            placed = [part.transform(frames[owners[key]]) if owners[key] in frames else part
                      for key, part in hardware.items()]
            hits = find_hits(placed, keys, pairs)
            if hits:
                findings.append(dict(frame=index, time_s=float(row["time_s"]), hits=hits))
            elif not findings:
                elbow = float(row["elbow_deg"])
                clear_peak = elbow if clear_peak is None else max(clear_peak, elbow)
            if index % 200 == 0:
                print("Rigid spring-rate check", index, "findings", len(findings), flush=True)
        filename = (f"nominal_spring_spool_{args.spool_rigid}_rigid_audit.json" if args.spool_rigid
                    else "return_spring_nominal_rate_rigid_audit.json")
        save(folder / filename, dict(
            generator=provenance(ROOT / "tools/audit_s1_series_elbow_spring_integration.py"),
            spring_rate_study=rate_ref, spool_study=spool_ref, trace=trace_ref, retention=retention,
            checked_pairs=len(pairs), sampled_frames=len(trace), clear_peak_elbow_deg=clear_peak,
            findings=findings, scope="All cross-body rigid pairs involving a moving body, including "
            "revised spring retention hardware. The deforming coil, steel cables, and TPU are excluded.",
        ))
        print("Checked", len(trace), "frames; findings", len(findings), "clear peak", clear_peak, flush=True)
        return
    added = {key: hardware[key] - original[key] if key in original else hardware[key] for key in revisions}
    fixed = bounded({key: part for key, part in hardware.items() if owners[key] not in MOBILE})
    moving = {key: part for key, part in hardware.items() if owners[key] in MOBILE}
    mf = A.I._manifold()
    start = np.array(mount["fixed_eye_source_mm"])
    end = np.array(mount["moving_eye_neutral_source_mm"])
    length = np.linalg.norm(end-start)
    axis = (end-start) / length
    pin_axis = np.array(mount["pin_axis_source"])
    spring_frame = np.c_[pin_axis, np.cross(axis, pin_axis), axis, start]
    spec = spring["design"]
    radius = .5 * spec["wire_diameter_mm"]
    mean = .5 * (spec["coil_OD_mm"] - spec["wire_diameter_mm"])
    turns = spec["active_turns"] + 2
    elbow = load_reference(load_reference(prior["passive"])["assembly"])
    muscles = {}
    for label, muscle in elbow["muscles"].items():
        definition = load_reference(muscle["source_definition"])
        path = np.array(muscle["wire_path_source_mm"])[1:-1]
        axis_m = np.array(definition["basis_columns_display"])[:, 2]
        gradient = axis_m / ((path[-1]-path[0]) @ axis_m)
        _, vertices, faces = indexed_solid(prior["mesh_sources"][label])
        muscles[label] = (muscle["anatomical_id"], path, weigh(path, path[0], gradient),
                          vertices, weigh(vertices, path[0], gradient), faces)
    findings = []
    for index, (row, checked) in enumerate(zip(trace, prior["records"], strict=True)):
        data.qpos[:] = [float(row[f"qpos_{i}"]) for i in range(model.nq)]
        mujoco.mj_kinematics(model, data)
        frames = transforms(model, data, neutral)
        h = frames["humerus"]
        travel = max(0., float(row["slide_mm"]))
        all_parts = {**fixed, **bounded({key: part.transform(frames[owners[key]]) for key, part in moving.items()})}
        changed = {key: all_parts[key] for key in revisions}
        # The baseline establishes the other rigid pairs. Use the complete
        # added volume for cross-body checks and also screen all new screws
        # against the parts on their own body, including their bores and pins.
        rigid_hits = {}
        for key, part in added.items():
            if part.is_empty() or part.volume() <= .001:
                continue
            placed = part.transform(frames[owners[key]])
            targets = {k: v for k, v in all_parts.items() if k != key
                       and (owners[k] != owners[key] or key.endswith("dog_screw"))}
            hit = overlaps(placed, targets)
            if hit:
                rigid_hits[key] = hit
        coil = sweep(spring_path(length+travel, turns, radius, mean), radius).transform(spring_frame).transform(h)
        coil_hits = overlaps(coil, all_parts)
        wire_targets = {**changed, **bounded({"series:return_spring": coil})}
        wires = {}
        tpu = {}
        for label, (head, path, weights, vertices, vertex_weights, faces) in muscles.items():
            distal = frames[f"{head}_distal_pivot"]
            r = h[:, :3].T
            t = h[:, 3]
            shift = r @ (frames["FJ1477_carriage"][:, 3] - t) if label == "extensor" else np.zeros(3)
            affine = r @ distal[:, :3], r @ (distal[:, 3] - t), shift
            points = displacement_map(path, weights, *affine) @ r + t
            hits = overlaps(A.D.outer_envelope(points, .30), wire_targets)
            if hits:
                wires[label] = hits
            if checked["TPU_sample"]:
                points = displacement_map(vertices, vertex_weights, *affine) @ r + t
                mesh = mf.Manifold(mf.Mesh64(np.ascontiguousarray(points), faces))
                if mesh.status() != mf.Error.NoError:
                    raise ValueError(f"Invalid interpolated {label} at frame {index}")
                hits = overlaps(mesh, wire_targets)
                if hits:
                    tpu[label] = hits
        if rigid_hits or coil_hits or wires or tpu:
            findings.append(dict(frame=index, time_s=float(row["time_s"]), rigid_hits_mm3=rigid_hits,
                                 spring_hits_mm3=coil_hits, wire_hits_mm3=wires, TPU_hits_mm3=tpu))
        if index % 200 == 0:
            print("Retained spring check", index, "findings", len(findings), flush=True)
    if args.spool_flexible:
        output = f"nominal_spring_spool_{spool_case}_flexible_audit.json"
    elif args.nominal_rate_flexible:
        output = "return_spring_nominal_rate_flexible_audit.json"
    elif args.anterior_friction:
        output = "return_spring_anterior_friction_integration.json"
    else:
        output = "return_spring_shoulder_elbow_integration.json"
    save(folder / output, dict(
        generator=provenance(ROOT / "tools/audit_s1_series_elbow_spring_integration.py"),
        retention=retention, model=extended["model"], trace=trace_ref, friction=friction_ref,
        spring_rate_study=rate_ref, spool_study=spool_ref, spool_rigid_audit=spool_rigid_ref,
        changed_sources=revisions,
        changed_owners={key: owners[key] for key in revisions}, nominal_rigid_and_wire_frames=len(trace),
        elbow_TPU_frames=sum(row["TPU_sample"] for row in prior["records"]), findings=findings,
        scope="Incremental geometry replay using a recorded dynamics trajectory."
              " Solid coil and revised rigid geometry checked; both elbow TPU shapes use prescribed cap interpolation."
              " Cable envelopes follow the interpolated nominal paths."
              " Load response, friction, pin or screw strength and shoulder TPU geometry are outside this check."))
    print("Checked", len(trace), "frames; findings", len(findings), flush=True)


if __name__ == "__main__":
    main()
