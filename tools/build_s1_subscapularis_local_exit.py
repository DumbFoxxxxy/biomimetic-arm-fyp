"""Screen a supported FJ1504 cap outlet with a curved local bore."""

import argparse
from pathlib import Path

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_elbow_series_guide import unit
from build_s1_teres_origin_mount import A, export, load_component, save
from study_elbow_series_package import bounded, overlaps


def outlet_path(outer, normal, lateral, radius, turn_deg, extension):
    turn = np.deg2rad(turn_deg)
    angles = np.linspace(0., turn, int(np.ceil(turn_deg / 5.)) + 1)
    arc = outer + radius * (np.outer(np.sin(angles), normal) + np.outer(1. - np.cos(angles), lateral))
    direction = np.cos(turn) * normal + np.sin(turn) * lateral
    return np.vstack([outer - .8 * normal, arc, arc[-1] + extension * direction]), direction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bidirectional", action="store_true")
    args = parser.parse_args()
    folder = DEST.parent / "integration_review/s1_series_elbow"
    assembly_file = folder / "assembly_scapular_relief.json"
    revision_file = folder / "scapular_thorax_relief.json"
    assembly, revision = map(read, (assembly_file, revision_file))
    if assembly["scapular_mount_revision"] != provenance(revision_file):
        raise RuntimeError("Scapular revision changed")
    route_file = ROOT / revision["route_reference"]["file"]
    if provenance(route_file) != revision["route_reference"]:
        raise RuntimeError("Origin route reference changed")
    route = next(row for row in read(route_file)["routes"] if row["head"] == "FJ1504")
    port = route["ports"]["origin"]
    shift = next(row["selected"]["shift_mm"] for row in revision["reports"] if row["head"] == "FJ1504")
    normal = unit(np.asarray(port["normal_source"]))
    outer = np.asarray(port["outer_face_source_mm"]) + shift * normal
    inner = np.asarray(port["inner_face_source_mm"]) + shift * normal
    first = unit(np.cross(normal, np.eye(3)[np.argmin(abs(normal))]))
    second = np.cross(normal, first)
    cap_name = "cap:FJ1504/origin"
    cap = load_component(assembly["sources"][cap_name])
    if args.bidirectional:
        assembly_file = folder / "assembly_FJ1500_front_mount.json"
        assembly = read(assembly_file)
    solids = {name: load_component(row) for name, row in assembly["sources"].items() if name != cap_name}
    obstacles = bounded(solids)
    thorax = bounded({key: solid for key, solid in solids.items() if key.startswith("thorax_reference:")})
    scapula = np.asarray(read(DEST / "elbow_passive_model/model_check.json")["neutral_origins_m"]["scapula"]) * 1000.
    clavicle = np.asarray(read(DEST / "elbow_passive_model/model_check.json")["neutral_origins_m"]["clavicle"]) * 1000.
    angles = (-2., -1.5, -1., -.5, .5, 1., 1.5, 2.) if args.bidirectional else (1., 3., 5., 10.)
    poses = [(joint, index, angle, A.I.rotation_about(pivot, axis, angle))
             for joint, pivot in (("AC", scapula), ("SC", clavicle))
             for index, axis in enumerate(np.eye(3)) for angle in angles]
    trials = []
    choices = []
    for radius in (6., 5., 4., 3., 2.5):
        for turn in (90., 75., 60.):
            for phase in range(0, 360, 15):
                angle = np.deg2rad(phase)
                lateral = np.cos(angle) * first + np.sin(angle) * second
                path, direction = outlet_path(outer, normal, lateral, radius, turn, 4.)
                envelope = A.D.outer_envelope(path, 1.65)
                row = dict(radius_mm=radius, turn_deg=turn, azimuth_deg=phase)
                trials.append(row)
                hits = overlaps(envelope - cap, obstacles)
                if hits:
                    row["neutral_overlap_mm3"] = hits
                    continue
                lumen = A.D.outer_envelope(np.vstack([inner - 1. * normal, path[1:], path[-1] + direction]), .86)
                outlet = ((cap + envelope) - lumen).as_original().simplify(1e-6)
                wire = A.D.outer_envelope(np.vstack([inner + .1 * normal, path[1:]]), .3)
                amount = float((wire ^ outlet).volume())
                if amount > .001 or sum(piece.volume() > 1e-8 for piece in outlet.decompose()) != 1:
                    row["disconnected_or_lumen_overlap_mm3"] = amount
                    continue
                assessed_part = outlet if args.bidirectional else outlet - cap
                sample_hits = [overlaps(assessed_part.transform(matrix), thorax) for *_, matrix in poses]
                row["shoulder_pose_hits"] = sum(bool(hit) for hit in sample_hits)
                row["first_pose_hit_abs_deg"] = min((abs(pose[2]) for pose, hit in zip(poses, sample_hits) if hit),
                                                    default=max(map(abs, angles)) + 1.)
                if not args.bidirectional:
                    row["positive_SC_AC_pose_hits"] = row["shoulder_pose_hits"]
                row["maximum_pose_overlap_mm3"] = max((sum(hit.values()) for hit in sample_hits), default=0.)
                choices.append((row, outlet, path, direction, wire))
    selected = None
    ranked = sorted(choices, key=lambda item: (
        -item[0]["first_pose_hit_abs_deg"] if args.bidirectional else 0.,
        item[0]["shoulder_pose_hits"], -item[0]["radius_mm"], item[0]["maximum_pose_overlap_mm3"]))
    variant = "_bidirectional" if args.bidirectional else ""
    for row, outlet, path, direction, wire in ranked:
        try:
            solid, source = export(outlet, folder / f"FJ1504_cap_curved_local_exit{variant}.obj", inner)
        except (AssertionError, ValueError) as error:
            row["export_error"] = str(error) or type(error).__name__
            continue
        clearance = overlaps(solid, obstacles)
        if clearance or float((wire ^ solid).volume()) > .001:
            row["readback_hits_mm3"] = clearance
            continue
        selected = dict(**row, source=source, path_source_mm=path.tolist(),
                        local_wire_path_source_mm=np.vstack([inner + .1 * normal, path[1:]]).tolist(),
                        exit_tangent_source=direction.tolist(), cap_inner_source_mm=inner.tolist(),
                        readback_wire_overlap_mm3=float((wire ^ solid).volume()))
        if args.bidirectional:
            selected["whole_cap_SC_AC_probes"] = [
                dict(joint=joint, source_axis_index=index, angle_deg=angle,
                     overlap_mm3=overlaps(solid.transform(matrix), thorax))
                for joint, index, angle, matrix in poses
            ]
        break
    save(folder / f"FJ1504_curved_local_exit{variant}.json", dict(
        generator=provenance(Path(__file__).resolve()), assembly=provenance(assembly_file),
        source_routes=provenance(route_file), outer_radius_mm=1.65, bore_radius_mm=.86,
        selected=selected, trials=trials, shoulder_pose_count=len(poses),
        scope="Cap-integrated curved local guide from the actual axial cap bore to an open end. "
        "Bidirectional probes use the complete outlet cap. Original probes use only added guide material. "
        "Neither cable bending endurance, drive, external connection, printability nor the TPU-side "
        "transition is established. Preserves previous rigid assembly as baseline."))
    print("lanes", len(trials), "neutral candidates", len(choices), "selected", selected, flush=True)


if __name__ == "__main__":
    main()
