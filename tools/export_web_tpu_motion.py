"""Pack the checked K1 and K7 TPU vertex motion for the interactive replay."""

import hashlib
import json

import mujoco
import numpy as np

from arm_action_timing import geometric_replay_indices
from audit_elbow_slider_short_guides import weigh
from audit_s1_series_elbow_integration import transforms
from audit_whole_arm_k7_tpu import affine, main as k7_replay
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_web_model import load_parts
from compare_opensim_whole_arm import load_k1
from refine_elbow_short_guide_relief import indexed_solid
from study_s1_local_force_geometry import load_reference


FOLDER = DEST.parent / "whole_arm_first_pass"
VIEWER = ROOT / "docs/viewer"
REPORT = FOLDER / "action_showcase/current_actions.json"
FIT_TOLERANCE_MM = .012
QUANTIZED_TOLERANCE_MM = .018


class K7Shapes:
    def __init__(self, folder, combo, sources, cuff_ref, cuffs, spring):
        audit = read(folder / "whole_arm_k7_tpu_moving_cuff_audit.json")
        if (audit["combined"] != combo or audit["tpu_sources"] != sources
                or audit["moving_cuff_reference"] != cuff_ref or spring):
            raise ValueError("K7 TPU inputs differ from the selected geometry audit")
        self.audit = audit["samples"]
        self.cuffs = cuffs
        self.meshes = {}
        self.wires = []
        self.angles = []
        self.poses = []
        self.routes = load_reference(load_reference(combo)["routed"])["routes"]

    def render(self, index, pose, muscles, paths):
        checked = self.audit[index]
        if checked["input_deg"] != pose["input_deg"] or checked["invalid_meshes"]:
            raise ValueError(f"Invalid K7 TPU sample {index}")
        for key, (vertices, faces) in muscles.items():
            if key not in self.meshes:
                self.meshes[key] = dict(faces=faces, frames=[])
            self.meshes[key]["frames"].append(vertices.astype(np.float32))
        self.angles.append(pose["humerus_direction_change_deg"])
        self.poses.append(pose["qpos"])
        self.wires.append(paths)

    def finish(self):
        if len(self.angles) != len(self.audit):
            raise ValueError("K7 capture missed part of the checked trajectory")


def k1_shapes(record, routes, sources):
    cap_report, model = load_k1()
    if cap_report["cap_relief"]["selected_source"] != record["cap"]:
        raise ValueError("K1 cap differs from the interactive rigid replay")
    data = mujoco.MjData(model)
    origins = {key: np.asarray(value) for key, value in
               read(FOLDER / "whole_arm_shared_scene_check.json")["neutral_body_origins_mm"].items()}
    action = record["actions"]["K1"]["playback"]
    outbound = action[:len(action) // 2 + 1]
    selected = sorted(set(geometric_replay_indices(len(action))))
    meshes, wires = {}, [{} for _ in selected]
    for role, wire in (("flexor", "FJ1486"), ("extensor", "FJ1477")):
        key = f"tpu:elbow_{role}"
        _, vertices, faces = indexed_solid(sources[key])
        route = routes[wire]
        points = np.asarray(route["points_source_mm"])
        first, last = points[[0, -1]]
        span = last - first
        weights = weigh(vertices, first, span / (span @ span))[:, None]
        route_weights = np.asarray(route["distal_weights"])[:, None]
        posed = []
        for index, source_index in enumerate(selected):
            data.qpos[:] = outbound[source_index]["qpos"]
            mujoco.mj_forward(model, data)
            frames = transforms(model, data, origins)
            a, b = (frames[owner] for owner in route["owners"])
            posed.append(((1. - weights) * affine(vertices, a) + weights * affine(vertices, b)).astype(np.float32))
            wires[index][wire] = (1. - route_weights) * affine(points, a) + route_weights * affine(points, b)
        meshes[key] = (faces, posed)
    return selected, meshes, wires


def fit_shapes(frames, base):
    # SVD of the small temporal Gram matrix. Store spatial modes only once;
    # the browser interpolates the corresponding scalar coefficients.
    matrix = np.asarray(frames, dtype=np.float64).reshape(len(frames), -1)
    matrix -= base.ravel()
    values, vectors = np.linalg.eigh(matrix @ matrix.T)
    vectors = vectors[:, np.argsort(values)[::-1]]
    residual = matrix.copy()
    modes = []
    for mode in range(len(frames)):
        if float(np.max(np.abs(residual))) <= FIT_TOLERANCE_MM:
            break
        temporal = vectors[:, mode]
        spatial = temporal @ matrix
        residual -= temporal[:, None] * spatial
        scale = float(np.max(np.abs(spatial))) / 32767.
        if scale == 0.:
            continue
        modes.append((np.rint(spatial / scale).astype("<i2"), scale, temporal))
    error = matrix.copy()
    for packed, scale, temporal in modes:
        error -= temporal[:, None] * (packed.astype(np.float32) * scale)
    maximum = float(np.max(np.abs(error)))
    if not np.isfinite(maximum) or maximum > QUANTIZED_TOLERANCE_MM:
        raise ValueError(f"Compressed TPU keyframes exceed error tolerance: {maximum:.5f} mm")
    return modes, maximum


def pack_array(buffer, data):
    offset = len(buffer)
    buffer.extend(data.tobytes())
    buffer.extend(b"\0" * (-len(buffer) % 4))
    return offset


def pack_mesh(buffer, vertices, faces, pivot):
    indices = faces.astype("<u2" if len(vertices) <= 65535 else "<u4")
    return dict(vertices=len(vertices), triangles=len(faces), index_bytes=indices.dtype.itemsize,
                position_offset=pack_array(buffer, ((vertices - pivot) * .001).astype("<f4")),
                index_offset=pack_array(buffer, indices))


def pack_cuffs(buffer, cuffs, pivot):
    hardware = {}
    for key, vertices, faces in load_parts(cuffs["moving_sources"]):
        poses = []
        for frame in cuffs["frames"]:
            basis = np.asarray(frame["matrices"][key])
            rotation = np.ascontiguousarray(basis[:, :3])
            np.testing.assert_allclose(rotation @ rotation.T, np.eye(3), atol=1e-8, rtol=0.)
            quaternion = np.empty(4)
            mujoco.mju_mat2Quat(quaternion, rotation.ravel())
            translation = (rotation @ pivot + basis[:, 3] - pivot) * .001
            poses.append([*translation, *quaternion[[1, 2, 3, 0]]])
        hardware[key] = dict(**pack_mesh(buffer, vertices, faces, pivot), poses=np.round(poses, 7).tolist())
    return hardware


def main():
    captures = []

    def capture(*args):
        replay = K7Shapes(*args)
        captures.append(replay)
        return replay

    k7_replay(["--render-replay", "--check-moving-cuffs"], replay_factory=capture)
    (k7,) = captures
    record = read(REPORT)
    playback = record["actions"]["K7"]["playback"]
    outward = playback[:len(playback) // 2 + 1]
    np.testing.assert_allclose([pose["qpos"] for pose in outward], k7.poses, atol=1e-10, rtol=0.)
    np.testing.assert_allclose([pose["input_value"] for pose in outward], k7.angles, atol=1e-10, rtol=0.)
    sources = read(FOLDER / "whole_arm_k7_tpu_moving_cuff_audit.json")["tpu_sources"]
    k1_keyframes, k1, k1_wires = k1_shapes(record, k7.routes, sources)
    rigid = read(VIEWER / "whole-arm-actions.json")
    if rigid["evidence"]["replay"] != provenance(REPORT):
        raise ValueError("Web rigid scene uses a different motion record")
    pivot = np.asarray(rigid["source_origin_mm"])
    parts, buffers = {}, {"elbow": bytearray(), "shoulder": bytearray()}
    for key, mesh in k7.meshes.items():
        elbow = key in k1
        chunk = "elbow" if elbow else "shoulder"
        buffer = buffers[chunk]
        base = mesh["frames"][0]
        faces = mesh["faces"]
        frames = mesh["frames"]
        if elbow:
            k1_faces, k1_frames = k1[key]
            np.testing.assert_array_equal(faces, k1_faces)
            frames = [*frames, *k1_frames]
        modes, error = fit_shapes(frames, base)
        geometry = pack_mesh(buffer, base, faces, pivot)
        packed_modes = [dict(offset=pack_array(buffer, encoded), scale_m=scale * .001)
                        for encoded, scale, _ in modes]
        coefficients = np.column_stack([time for _, _, time in modes]).tolist() if modes else [[] for _ in frames]
        parts[key] = dict(chunk=chunk, **geometry, modes=packed_modes,
                          coefficients={"K7": coefficients[:len(k7.angles)],
                                        **({"K1": coefficients[len(k7.angles):]} if elbow else {})},
                          maximum_coordinate_error_mm=error)
        print(key, len(modes), f"keyframe error {error:.4f} mm", flush=True)

    def wire_frames(frames):
        return [{key: np.round((points - pivot) * .001, 7).tolist() for key, points in frame.items()}
                for frame in frames]

    web_motion = read(VIEWER / rigid["motion_file"])
    actions = {}
    for key, keyframes, wires in (("K1", k1_keyframes, k1_wires),
                                  ("K7", list(range(len(k7.angles))), k7.wires)):
        raw = record["actions"][key]["playback"]
        indices = geometric_replay_indices(len(raw))
        np.testing.assert_allclose([sample[0] for sample in web_motion["actions"][key]["samples"]],
                                   [raw[index]["input_value"] for index in indices], atol=5.1e-6, rtol=0.)
        actions[key] = dict(keyframes=keyframes, indices=indices, wires=wire_frames(wires))
    np.testing.assert_allclose([frame["input_deg"] for frame in k7.cuffs["frames"]],
                               [frame["input_deg"] for frame in k7.audit], atol=1e-10, rtol=0.)
    hardware = pack_cuffs(buffers["shoulder"], k7.cuffs, pivot)
    chunks = {}
    for key, buffer in buffers.items():
        path = VIEWER / f"whole-arm-tpu-{key}.bin"
        path.write_bytes(buffer)
        chunks[key] = dict(file=path.name, bytes=len(buffer), sha256=hashlib.sha256(buffer).hexdigest())
    manifest = dict(mesh_sha256=rigid["mesh_sha256"], action_report=provenance(REPORT),
                    parts=parts, chunks=chunks, actions=actions, moving_cuffs=hardware)
    (VIEWER / "whole-arm-tpu.json").write_text(
        json.dumps(manifest, separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
