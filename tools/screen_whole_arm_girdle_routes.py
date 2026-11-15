"""Screen sampled scapular muscle chords on the shared first-pass assembly."""

from itertools import product

import mujoco
import numpy as np

from build_channel_first_girdle_candidates import representative
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, rod, save
from project_runtime import load_mujoco_model
from study_elbow_series_package import bounded, overlaps


def sample_patch(atlas_row, role, structure, folder, divisions, index, direction):
    record = next(row for row in atlas_row["regions"]
                  if row["role"] == role and row["structure"] == structure)
    path = folder / record["data"]
    with np.load(path) as archive:
        patch = {key: archive[key] for key in ("bone_points_source_mm", "muscle_face_areas_mm2")}
    all_points = patch["bone_points_source_mm"]
    indices = np.array_split(np.argsort(all_points[:, 2]), divisions)[index]
    patch = {key: values[indices] for key, values in patch.items()}
    points = patch["bone_points_source_mm"]
    ordered = np.argsort(points @ direction, kind="stable")
    # Keep interior samples rather than the extreme vertices of a proximity patch.
    selected = [
        representative(patch),
        points[ordered[round(.20 * (len(points) - 1))]].tolist(),
        points[ordered[round(.80 * (len(points) - 1))]].tolist(),
    ]
    return np.unique(selected, axis=0, return_index=True), provenance(path)


class SharedAssembly:
    def __init__(self, folder):
        self.check_file = folder / "geometry_mass_body_chain_check.json"
        check = read(self.check_file)
        self.assembly_file = ROOT / check["reference"]["file"]
        self.model_file = ROOT / check["generated_model"]["file"]
        for path, expected in ((self.assembly_file, check["reference"]),
                               (self.model_file, check["generated_model"])):
            if provenance(path) != expected:
                raise ValueError(f"Shared assembly input changed: {path}")
        self.model = load_mujoco_model(self.model_file)
        self.data = mujoco.MjData(self.model)
        self.origins = {k: np.asarray(v) * 1000 for k, v in check["neutral_body_origins_m"].items()}
        self.origins["thorax"] = np.zeros(3)
        assembly = read(self.assembly_file)
        sources = {**assembly["sources"], **assembly["joint_context_sources"]}
        self.solids = {row["source_key"]: (row["owner"], load_component(sources[row["source_key"]]))
                       for row in check["geom_owners"]}
        self.sc = self.model.joint("SC").dofadr.item()
        self.ac = self.model.joint("AC").dofadr.item()

    def pose(self, axis=0, sc=0., ac=(0., 0., 0.)):
        mujoco.mj_resetData(self.model, self.data)
        delta = np.zeros(self.model.nv)
        delta[self.sc + axis] = np.deg2rad(sc)
        delta[self.ac:self.ac + 3] = np.deg2rad(ac)
        mujoco.mj_integratePos(self.model, self.data.qpos, delta, 1.)
        mujoco.mj_forward(self.model, self.data)
        frames = {}
        for name, origin in self.origins.items():
            body = self.data.body(name)
            rotation = body.xmat.reshape(3, 3).copy()
            frames[name] = (rotation, body.xpos * 1000 - rotation @ origin)
        solids = {key: part.transform(np.column_stack(frames[owner]))
                  for key, (owner, part) in self.solids.items()}
        return frames, bounded(solids), self.data.qpos.copy()


def screen(endpoints, owners, pose, measure_gap=False):
    frames, solids, _ = pose
    a, b = [frames[owner][0] @ point + frames[owner][1]
            for point, owner in zip(endpoints, owners)]
    wire = rod(a, b, .30)
    row = dict(length_mm=float(np.linalg.norm(b - a)), rigid_hits_mm3=overlaps(wire, solids))
    if measure_gap:
        # Endpoints are on bone surfaces. Measure clearance of the interior only.
        direction = (b - a) / row["length_mm"]
        span = rod(a + direction * 2., b - direction * 2., .30)
        box = np.asarray(span.bounding_box())
        row["interior_gap_capped_3mm"] = min((
            A.I.gap(span, solid, 3.) for solid, other in solids.values()
            if np.all(box[:3] - 3. <= other[3:]) and np.all(other[:3] - 3. <= box[3:])
        ), default=3.)
    return row


def stand_off_candidate(candidates, owners, neutral):
    """Try parallel lifts of a short list of anatomical chords."""
    closest = sorted(candidates, key=lambda row: sum(row["rigid_hits_mm3"].values()))
    attempts = []
    for base in closest[:3]:
        points = np.asarray(base["source_points_mm"])
        axis = points[1] - points[0]
        axis /= np.linalg.norm(axis)
        side = np.cross(axis, np.eye(3)[np.argmin(abs(axis))])
        side /= np.linalg.norm(side)
        up = np.cross(axis, side)
        for distance in (2., 4., 6.):
            for angle in np.arange(8) * np.pi / 4:
                offset = distance * (np.cos(angle) * side + np.sin(angle) * up)
                new_points = points + offset
                hit_checks = screen(new_points, owners, neutral)
                trial = dict(base_sample_indices=base["sample_indices"],
                             bone_attach_reference_mm=points.tolist(), offset_source_mm=offset.tolist(),
                             source_points_mm=new_points.tolist(), **hit_checks)
                attempts.append(trial)
                if not trial["rigid_hits_mm3"]:
                    # The proposed seat axis also needs a path from each base to its port.
                    mounts = [overlaps(rod(a, b, .75), neutral[1])
                              for a, b in zip(points, new_points)]
                    trial["seat_axis_probe_hits_mm3"] = mounts
                    unintended = [{key: volume for key, volume in hits.items()
                                   if solids_owner_key(key, owner) is False}
                                  for owner, hits in zip(owners, mounts)]
                    if not any(unintended):
                        trial["seat_probe_note"] = "Nominal 1.5 mm probes; same-owner overlap is listed for follow-up."
                        return trial, attempts
    return None, attempts


def solids_owner_key(key, owner):
    return key == f"bone:{owner}" or (owner == "thorax" and key.startswith("thorax_reference:"))


def main():
    folder = DEST / "unified_SC_AC"
    shared = SharedAssembly(folder)
    contact_file = folder / "contact_following_paths.json"
    contact = read(contact_file)
    if contact["reference"] != provenance(shared.assembly_file):
        raise ValueError("Shoulder-girdle poses use a different assembly")
    patch_folder = DEST.parent / "anatomical_concept/attachment_regions"
    atlas_file = patch_folder / "atlas_and_regions.json"
    atlas = {row["id"]: row for row in read(atlas_file)["muscles"]}
    input_file = DEST / "upper_girdle_endpoint_candidates.json"
    inputs = read(input_file)["routes"]
    divisions = {head: sum(item["anatomical_id"] == head for item in inputs) for head in atlas}
    neutral = shared.pose()
    routes = []
    for item in inputs:
        owners = (item["anatomical_origin_body"], item["distal_body"])
        origin_structure = item["candidate_origin_structure"]
        if origin_structure == "cervical_spine":
            routes.append(dict(id=item["id"], muscle=item["module"],
                               status="missing_cervical_support",
                               anatomical_reference_source_mm=item["candidate_endpoints_source_mm"]))
            continue
        head = item["anatomical_id"]
        index = int(item["id"].rsplit("g", 1)[1]) - 1
        a0, b0 = item["candidate_endpoints_source_mm"]
        axis = np.asarray(b0) - a0
        axis /= np.linalg.norm(axis)
        sampled = [sample_patch(atlas[head], role, structure, patch_folder, divisions[head], index, axis)
                   for role, structure in (("origin", origin_structure), ("insertion", owners[1]))]
        points = [values[np.argsort(order)] for (values, order), _ in sampled]
        candidates = []
        for (a_index, a), (b_index, b) in product(enumerate(points[0]), enumerate(points[1])):
            row = screen((a, b), owners, neutral, measure_gap=True)
            candidates.append(dict(sample_indices=[a_index, b_index],
                                   source_points_mm=[a.tolist(), b.tolist()],
                                   body_local_points_mm=[(point - shared.origins[owner]).tolist()
                                                         for point, owner in zip((a, b), owners)],
                                   **row))
        selected = next((row for row in candidates if not row["rigid_hits_mm3"]), None)
        attempts = []
        if selected is None:
            selected, attempts = stand_off_candidate(candidates, owners, neutral)
        if selected:
            selected["body_local_points_mm"] = [(np.asarray(point) - shared.origins[owner]).tolist()
                                                for point, owner in zip(selected["source_points_mm"], owners)]
        routes.append(dict(id=item["id"], muscle=item["module"], owners=owners,
                           patches=[entry[1] for entry in sampled], candidates=candidates,
                           stand_off_trials=attempts,
                           neutral_selection=selected,
                           status="chord_candidate_pending_seats" if selected else "needs_local_guide_or_seat",
                           motion=[]))
        print(item["id"], sum(not row["rigid_hits_mm3"] for row in candidates), "clear anatomical chords;",
              "offset port candidate" if selected and attempts else routes[-1]["status"], flush=True)

    for branch in contact["paths"]:
        eligible = [row for row in routes if row.get("neutral_selection")]
        if not eligible:
            break
        axis = branch["axis"]
        states = [{"SC_angle_deg": 0., "AC_rotvec_deg": [0., 0., 0.]}, *branch["states"]]
        poses = []
        for route in eligible:
            route["motion"].append(dict(axis=axis, direction=branch["direction"], outgoing=[],
                                        return_samples=0, first_limit=branch["first_failure"]))
        for state in states:
            pose = shared.pose(axis, state["SC_angle_deg"], state["AC_rotvec_deg"])
            poses.append(pose)
            for route in eligible[:]:
                trial = screen(route["neutral_selection"]["source_points_mm"], route["owners"], pose,
                               measure_gap=True)
                result = route["motion"][-1]
                angle = state["SC_angle_deg"]
                if trial["rigid_hits_mm3"]:
                    result["first_limit"] = dict(reason="wire_interference", angle_deg=angle,
                                                 hits=trial["rigid_hits_mm3"])
                    eligible.remove(route)
                    continue
                result["outgoing"].append(dict(angle_deg=angle, ac_rotvec_deg=state["AC_rotvec_deg"], **trial))
        for route in routes:
            if not route.get("neutral_selection"):
                continue
            result = route["motion"][-1]
            rows = result["outgoing"]
            for row, pose in zip(reversed(rows), reversed(poses[:len(rows)]), strict=True):
                trial = screen(route["neutral_selection"]["source_points_mm"], route["owners"], pose)
                if trial["rigid_hits_mm3"] or not np.isclose(trial["length_mm"], row["length_mm"]):
                    raise ValueError(f"Return check differs for {route['id']}")
                result["return_samples"] += 1
            angles = np.deg2rad([row["angle_deg"] for row in rows])
            lengths = np.array([row["length_mm"] for row in rows])
            result["length_excursion_mm"] = float(np.ptp(lengths)) if len(lengths) else None
            # Travel per radian along the coordinated motion, not an isolated joint moment arm.
            result["effective_leverage_along_path_mm"] = (-np.diff(lengths) / np.diff(angles)).tolist()

    output = DEST.parent / "whole_arm_first_pass"
    output.mkdir(exist_ok=True)
    save(output / "girdle_route_screen.json", dict(
        assembly=provenance(shared.assembly_file), model=provenance(shared.model_file),
        coordinated_paths=provenance(contact_file), candidate_references=provenance(input_file),
        attachment_atlas=provenance(atlas_file), generator=provenance(ROOT / "tools/screen_whole_arm_girdle_routes.py"),
        rigid_solid_count=len(shared.solids), wire_diameter_mm=.60, interference_threshold_mm3=.001,
        scope="Flat-ended straight wire envelopes. Interior clearance excludes 2 mm at each end. "
              "Offset port and seat-axis probes are candidate geometry; actual reaction seats and TPU are pending. "
              "SC/AC trajectories are prescribed "
              "contact-following samples; signed leverage is reported along SC angle with coordinated AC, "
              "not a measured motion range or isolated moment arm.",
        routes=routes))


if __name__ == "__main__":
    main()
