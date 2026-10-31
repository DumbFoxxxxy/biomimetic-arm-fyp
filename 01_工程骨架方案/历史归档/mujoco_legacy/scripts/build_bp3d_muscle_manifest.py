"""Fuse the BodyParts3D muscles into drive channels and extract wire paths.

The printed X-lattice muscle wraps the drive wire, so the wire has to run
along the centre of the muscle volume.  For every channel the BodyParts3D
muscle meshes are sliced perpendicular to their principal axis and the
centroid of each cross-section becomes a point on the wire path.  A channel
with more wires than muscle elements gets extra fibres by cutting the largest
element in half along its second principal axis.

The output is ``config/bp3d_muscle_manifest.json``; every path records its two
attachment points, the full centre line for the CAD side, the tendon sites the
MuJoCo builder needs, and a report of how close the wire passes to each bone.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import bp3d_common as B


def bone_surfaces() -> tuple[dict[str, cKDTree], dict[str, dict[str, np.ndarray]]]:
    """Per-bone search structures for nearest-point and inside tests."""

    trees: dict[str, cKDTree] = {}
    surfaces: dict[str, dict[str, np.ndarray]] = {}
    for spec in B.BONE_SPECS:
        name = spec["name"]
        vertices, faces = B.read_obj(B.ASSET_DIR / f"{name}.obj")
        tri = vertices[faces]
        samples = tri.reshape(-1, 3)
        if len(samples) > 80000:
            samples = samples[:: int(np.ceil(len(samples) / 80000))]
        trees[name] = cKDTree(samples)
        normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        lengths = np.linalg.norm(normals, axis=1)
        keep = lengths > 1e-12
        surfaces[name] = {
            "centroid": tri[keep].mean(axis=1),
            "normal": normals[keep] / lengths[keep, None],
            "normal_tree": cKDTree(tri[keep].mean(axis=1)),
        }
    return trees, surfaces


def load_element_triangles(file_id: str) -> np.ndarray:
    vertices, faces = B.read_obj(B.element_path(file_id))
    return vertices[faces] * B.MM


def fibres_for_channel(
    channel: dict,
    rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    """One entry per drive wire: the muscle volume and the lateral offset.

    A channel with more wires than muscle elements gets extra parallel fibres
    inside the largest muscle.  The fibres share the centre line but are
    displaced sideways, which is how several Bowden cables sit side by side in
    one printed muscle belly.  Cutting the mesh in half instead would leave a
    fan-shaped muscle such as the serratus anterior with no scapular
    attachment on one of the halves.
    """

    elements: list[tuple[str, np.ndarray, float]] = []
    for pattern in channel["elements"]:
        for file_id in B.resolve_elements(rows, pattern):
            triangles = load_element_triangles(file_id)
            volume = float(B.mesh_mass_properties(triangles, 1.0)["volume_m3"])
            elements.append((file_id, triangles, volume))
    if not elements:
        raise FileNotFoundError(f"No muscle element found for {channel['id']}")
    elements.sort(key=lambda item: item[2], reverse=True)

    wanted = int(channel.get("channels", 1))
    fibres: list[dict[str, object]] = [
        {
            "element": file_id,
            "triangles": triangles,
            "part": "whole",
            "index": index,
            "offset_m": 0.0,
        }
        for index, (file_id, triangles, _) in enumerate(elements[:wanted])
    ]
    if len(fibres) < wanted:
        largest = max(fibres, key=lambda item: _volume(item["triangles"]))
        triangles = largest["triangles"]
        vertices = triangles.reshape(-1, 3)
        _, axes = B.principal_frame(vertices)
        lateral = axes[:, 1]
        width = float(np.ptp((vertices - vertices.mean(axis=0)) @ lateral))
        step = float(np.clip(width * 0.18, 0.002, 0.008))
        extra = wanted - len(fibres)
        for number in range(extra):
            sign = 1.0 if number % 2 == 0 else -1.0
            level = sign * (1.0 + number // 2)
            fibres.append(
                {
                    "element": largest["element"],
                    "triangles": triangles,
                    "part": f"fibre{number + 1}",
                    "index": largest["index"],
                    "offset_m": level * step,
                    "lateral": lateral,
                }
            )
    return fibres


def _volume(triangles: np.ndarray) -> float:
    return float(B.mesh_mass_properties(triangles, 1.0)["volume_m3"])


def _query(trees: dict[str, cKDTree], bone: str, point: np.ndarray) -> float:
    return float(trees[bone].query(point)[0])


def end_points(
    triangles: np.ndarray,
    trees: dict[str, cKDTree],
    origin_bone: str,
    insertion_bone: str,
    fallback_bones: tuple[str, ...],
    fraction: float = 0.14,
    max_residual_m: float = 0.035,
) -> dict[str, object]:
    """Attachment points of one muscle: which bone each end sits on, and where."""

    vertices = triangles.reshape(-1, 3)
    centre, axes = B.principal_frame(vertices)
    projection = (vertices - centre) @ axes[:, 0]
    low, high = np.percentile(projection, [fraction * 100.0, (1.0 - fraction) * 100.0])
    ends = {
        "a": vertices[projection <= low].mean(axis=0),
        "b": vertices[projection >= high].mean(axis=0),
    }

    def snap(bone: str, label: str) -> dict[str, object]:
        distance, index = trees[bone].query(ends[label])
        return {
            "bone": bone,
            "point_m": np.asarray(trees[bone].data[index], dtype=float),
            "distance_m": float(distance),
        }

    # Roles are decided end by end against the two declared bones: an end is
    # the origin when it is nearer the origin bone than the insertion bone.
    # The declared pair is kept because the subscapularis lies on the anterior
    # scapula only a few millimetres from the ribs, where a pure nearest-bone
    # rule would bolt it to the thorax.
    origin_is_a = _query(trees, origin_bone, ends["a"]) <= _query(
        trees, insertion_bone, ends["a"]
    )
    origin_is_b = _query(trees, origin_bone, ends["b"]) <= _query(
        trees, insertion_bone, ends["b"]
    )
    if origin_is_a and not origin_is_b:
        origin, insertion = snap(origin_bone, "a"), snap(insertion_bone, "b")
    elif origin_is_b and not origin_is_a:
        origin, insertion = snap(origin_bone, "b"), snap(insertion_bone, "a")
    else:
        straight = _query(trees, origin_bone, ends["a"]) + _query(
            trees, insertion_bone, ends["b"]
        )
        swapped = _query(trees, origin_bone, ends["b"]) + _query(
            trees, insertion_bone, ends["a"]
        )
        if straight <= swapped:
            origin, insertion = snap(origin_bone, "a"), snap(insertion_bone, "b")
        else:
            origin, insertion = snap(origin_bone, "b"), snap(insertion_bone, "a")

    if max(origin["distance_m"], insertion["distance_m"]) > max_residual_m:
        scored = []
        for label, point in ends.items():
            for bone in fallback_bones:
                distance, index = trees[bone].query(point)
                scored.append(
                    {
                        "label": label,
                        "bone": bone,
                        "distance_m": float(distance),
                        "point_m": np.asarray(trees[bone].data[index], dtype=float),
                    }
                )
        scored.sort(key=lambda item: item["distance_m"])
        best = scored[0]
        other = next(item for item in scored if item["label"] != best["label"])
        origin, insertion = best, other

    return {
        "origin": {
            "bone": origin["bone"],
            "point_m": origin["point_m"],
            "distance_m": origin["distance_m"],
        },
        "insertion": {
            "bone": insertion["bone"],
            "point_m": insertion["point_m"],
            "distance_m": insertion["distance_m"],
        },
    }


def resample(points: np.ndarray, count: int) -> np.ndarray:
    """Arc-length resample a polyline, keeping both endpoints."""

    if len(points) <= count:
        return points
    segments = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cumulative = np.concatenate([[0.0], np.cumsum(segments)])
    targets = np.linspace(0.0, cumulative[-1], count)
    out = [points[0]]
    for target in targets[1:-1]:
        index = int(np.clip(np.searchsorted(cumulative, target) - 1, 0, len(points) - 2))
        span = max(cumulative[index + 1] - cumulative[index], 1e-12)
        weight = (target - cumulative[index]) / span
        out.append(points[index] + weight * (points[index + 1] - points[index]))
    out.append(points[-1])
    return np.asarray(out)


def wire_path(
    triangles: np.ndarray,
    attachments: dict[str, object],
    trees: dict[str, cKDTree],
    sites: int = 8,
    offset_m: float = 0.0,
    lateral: np.ndarray | None = None,
) -> dict[str, object]:
    """Centre line from the origin attachment to the insertion attachment.

    A non-zero ``offset_m`` shifts the interior of the line sideways, so a
    second wire can share the same muscle belly without leaving it.
    """

    line = B.centre_line(triangles, samples=24)
    origin = np.asarray(attachments["origin"]["point_m"], dtype=float)
    insertion = np.asarray(attachments["insertion"]["point_m"], dtype=float)
    if np.linalg.norm(line[0] - origin) > np.linalg.norm(line[-1] - origin):
        line = line[::-1]
    full = np.vstack([origin, line, insertion])
    # Drop centre-line points that sit behind an attachment.
    keep = np.linalg.norm(np.diff(full, axis=0), axis=1) > 1e-6
    full = np.vstack([full[0], full[1:][keep]])
    interior = resample(full, sites)
    interior[0] = origin
    interior[-1] = insertion
    if offset_m and lateral is not None and len(interior) > 2:
        interior[1:-1] = interior[1:-1] + offset_m * np.asarray(lateral, dtype=float)

    origin_bone = str(attachments["origin"]["bone"])
    insertion_bone = str(attachments["insertion"]["bone"])
    bones = [origin_bone, insertion_bone]
    distances = np.vstack([trees[bone].query(interior)[0] for bone in bones])
    owners = np.argmin(distances, axis=0)
    # The path is monotone: a point may not flip back to the origin bone after
    # it has moved onto the insertion bone.
    changes = np.flatnonzero(np.diff(owners) == 1)
    if len(changes):
        owners[int(changes[0]) + 1 :] = 1
    owners[0] = 0
    owners[-1] = 1

    sites_payload = [
        {
            "body": bones[int(owner)],
            "pos_m": point,
            "role": "origin" if index == 0 else ("insertion" if index == len(interior) - 1 else "via"),
        }
        for index, (point, owner) in enumerate(zip(interior, owners))
    ]
    length = float(np.linalg.norm(np.diff(interior, axis=0), axis=1).sum())
    return {
        "polyline_m": full,
        "sites": sites_payload,
        "length_m": length,
        "origin_bone": origin_bone,
        "insertion_bone": insertion_bone,
    }


def clearance_report(
    polyline: np.ndarray,
    surfaces: dict[str, dict[str, np.ndarray]],
    bodies: list[str],
    step: float = 0.002,
    far_mm: float = 80.0,
) -> dict[str, object]:
    """How close the wire passes to each bone, and whether it enters one.

    Bones further away than ``far_mm`` from the whole path are reported as
    clear without a full test, which keeps the thorax from dominating the run
    time for a wire that only crosses the forearm.
    """

    segments = np.linalg.norm(np.diff(polyline, axis=0), axis=1)
    samples = [polyline[0]]
    for index, length in enumerate(segments):
        count = max(1, int(np.ceil(length / step)))
        for k in range(1, count + 1):
            weight = k / count
            samples.append(
                polyline[index] + weight * (polyline[index + 1] - polyline[index])
            )
    samples = np.asarray(samples)

    report: dict[str, object] = {}
    for name, entry in surfaces.items():
        centroids = entry["centroid"]
        distance, index = entry["normal_tree"].query(samples, workers=-1)
        nearest = centroids[index]
        signed = np.einsum("ij,ij->i", samples - nearest, entry["normal"][index])
        if float(distance.min()) > far_mm * B.MM:
            report[name] = {
                "min_distance_mm": round(float(distance.min()) * 1000.0, 2),
                "samples_inside": 0,
                "samples_near": 0,
            }
            continue
        # Clearly behind the nearest triangle means the sample is inside the
        # bone volume.  The band keeps a grazing wire from being flagged.
        inside = (signed < -0.004) & (distance < 0.02)
        report[name] = {
            "min_distance_mm": round(float(distance.min()) * 1000.0, 2),
            "samples_inside": int(np.count_nonzero(inside)),
            "samples_near": int(np.count_nonzero(distance <= far_mm * B.MM)),
        }
    others = [
        entry["min_distance_mm"]
        for name, entry in report.items()
        if name not in bodies
    ]
    report["path_bones"] = bodies
    report["worst_other_bone_mm"] = round(float(min(others)) if others else float("nan"), 2)
    report["worst_own_bone_mm"] = round(
        float(min(report[bone]["min_distance_mm"] for bone in bodies if bone in report)),
        2,
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=B.MUSCLE_MANIFEST_PATH)
    parser.add_argument("--sites", type=int, default=8)
    parser.add_argument("--include-optional", action="store_true")
    parser.add_argument("--skip-clearance", action="store_true")
    args = parser.parse_args()

    rows = B.load_index()
    trees, surfaces = bone_surfaces()
    all_bones = tuple(spec["name"] for spec in B.BONE_SPECS)
    muscle_dir = B.ASSET_DIR / "muscles"
    muscle_dir.mkdir(parents=True, exist_ok=True)

    contract = B.load_json(B.CONTRACT_PATH)
    channels = list(B.CHANNELS)
    if args.include_optional:
        channels += list(B.OPTIONAL_CHANNELS)

    manifest_channels: list[dict[str, object]] = []
    print(f"building {len(channels)} drive channels")
    for channel in channels:
        preferred = (channel["origin_bone"], channel["insertion_bone"])
        candidates = preferred + tuple(
            name for name in all_bones if name not in preferred
        )
        fibres = fibres_for_channel(channel, rows)
        paths: list[dict[str, object]] = []
        for index, fibre in enumerate(fibres):
            triangles = np.asarray(fibre["triangles"], dtype=float)
            attachments = end_points(
                triangles, trees, preferred[0], preferred[1], candidates
            )
            path = wire_path(
                triangles,
                attachments,
                trees,
                sites=args.sites,
                offset_m=float(fibre.get("offset_m", 0.0)),
                lateral=fibre.get("lateral"),
            )
            body_pair = [path["origin_bone"], path["insertion_bone"]]
            if not args.skip_clearance:
                path["clearance"] = clearance_report(
                    np.asarray(path["polyline_m"]), surfaces, body_pair
                )
            path["id"] = f"{channel['id']}_{index}"
            path["element"] = fibre["element"]
            path["part"] = fibre["part"]
            path["volume_cm3"] = round(_volume(triangles) * 1e6, 1)
            path["origin_distance_mm"] = round(
                float(attachments["origin"]["distance_m"]) * 1000.0, 2
            )
            path["insertion_distance_mm"] = round(
                float(attachments["insertion"]["distance_m"]) * 1000.0, 2
            )
            path["expected_origin_bone"] = preferred[0]
            path["expected_insertion_bone"] = preferred[1]
            B.write_obj(
                muscle_dir / f"{channel['id']}_{index}.obj",
                triangles.reshape(-1, 3),
                _unique_faces(triangles),
                header=f"{channel['id']} wire {index} muscle volume",
            )
            paths.append(path)
            worst = path.get("clearance", {}).get("worst_other_bone_mm", float("nan"))
            inside = sum(
                entry.get("samples_inside", 0)
                for name, entry in path.get("clearance", {}).items()
                if isinstance(entry, dict)
            )
            print(
                f"  {path['id']:28s} {path['origin_bone']:9s} -> "
                f"{path['insertion_bone']:9s} len={path['length_m']*1000:6.1f} mm "
                f"vol={path['volume_cm3']:7.1f} cm3  worst other bone={worst:6.2f} mm"
                f"  inside={inside}"
            )
            if path["origin_bone"] == path["insertion_bone"]:
                print(f"    ! {path['id']} attaches to the same bone twice")
        manifest_channels.append(
            _finalise_channel(channel, paths, preferred, contract)
        )

    manifest = {
        "schema_version": 1,
        "model_name": "right_arm_bp3d_cable_drive_v1",
        "units": "meters",
        "source": str(B.BP3D_DIR),
        "defaults": {
            "timestep": 0.001,
            "gravity": [0.0, 0.0, -9.81],
            "integrator": "implicitfast",
            "slider_offset_m": 0.09,
            "slider_radius_m": 0.008,
            "slider_kp": 20000.0,
            "slider_kv": 400.0,
            "release_m": 0.06,
            "wire_width_m": 0.0015,
            "pretension_n": 0.2,
            "reference_slack_m": 0.004,
            "damping_ratio": 0.7,
        },
        "notes": (
            "Wire paths follow the centre line of the BodyParts3D muscle "
            "volumes because the printed TPU X-lattice wraps the wire. "
            "lattice_muscle / lattice_tendon are the Yoshimura 2026 cell "
            "parameters for the deformable belly and the stiff ends."
        ),
        "lattice_reference": {
            "muscle_compressive_stiffness_n_per_mm": 1.3,
            "tendon_tensile_stiffness_n_per_mm": 15.0,
            "source": "Yoshimura, Kawaharazuka, Okada 2026, Figure 4",
        },
        "channels": manifest_channels,
    }
    B.save_json(args.manifest, manifest)
    total_paths = sum(entry["channels"] for entry in manifest_channels)
    print(f"\n{len(manifest_channels)} channels, {total_paths} wires")
    dropped_total = sum(len(entry["dropped_paths"]) for entry in manifest_channels)
    if dropped_total:
        print(f"dropped {dropped_total} mis-routed wires:")
        for entry in manifest_channels:
            for item in entry["dropped_paths"]:
                print(
                    f"  {item['id']:30s} dL/dq = {item['moment_arm_mm_per_rad']:8.2f} "
                    f"mm/rad on {entry['id']}"
                )
    conflicts = [entry["id"] for entry in manifest_channels if entry["sign_conflict"]]
    if conflicts:
        print(f"channels with a sign conflict: {conflicts}")
    print(f"manifest: {args.manifest}")


def _series_stiffness(paths: list[dict[str, object]]) -> float:
    """Series elasticity of one wire, from the lattice cell measurements.

    Each deformable muscle cell behaves as a 1.3 N/mm spring and the stiff
    end cells as 15 N/mm springs.  The belly and the two tendon-like ends are
    mechanically in series, so the softest element dominates:
    ``1/k = 1/k_belly + 1/k_end``.
    """

    if not paths:
        return 0.0
    length = float(np.mean([path["length_m"] for path in paths]))
    cells = max(1.0, length / (B.LATTICE_MUSCLE["x_mm"] * B.MM))
    belly = 1.3e3 / cells
    ends = 15.0e3 / max(1.0, cells / 4.0)
    return float(1.0 / (1.0 / belly + 1.0 / ends))


def _finalise_channel(
    channel: dict,
    paths: list[dict[str, object]],
    preferred: tuple[str, str],
    contract: dict,
) -> dict[str, object]:
    """Keep only the wires that pull the channel's joint the intended way.

    A muscle group is one printed body wrapping all of its wires, so every
    wire in the group has to pull the same way.  A wide fan such as the
    serratus anterior spans both sides of the joint axis, and the fibres on
    the wrong side would cancel the others, so they are dropped and reported.
    """

    target = B.CHANNEL_TARGETS.get(channel["id"])
    dropped: list[dict[str, object]] = []
    keep_mask: list[bool] = []
    if target is not None:
        joint, axis_name, action_sign = target
        wanted = -action_sign
        radius = float(contract["joints"][joint]["design_radius_m"]) + 0.003
        for path in paths:
            # Give the wire a wrapping point on the parent bone so its moment
            # arm cannot change sign inside the joint's range.
            sites, pulley = B.insert_joint_pulley(
                path["sites"], contract, joint, axis_name, radius, wanted
            )
            if pulley is not None:
                path["sites"] = sites
                path["pulley"] = pulley
        for path in paths:
            moment = B.path_moment_arm(path["sites"], contract, joint, axis_name)
            path["moment_arm_mm_per_rad"] = round(moment * 1000.0, 2)
            path["moment_arm_profile_mm_per_rad"] = [
                round(value * 1000.0, 2)
                for value in B.wire_moment_arm_profile(
                    path["sites"], contract, joint, axis_name, samples=7
                )
            ]
            if moment * wanted <= 0.0:
                dropped.append(
                    {
                        "id": path["id"],
                        "moment_arm_mm_per_rad": path["moment_arm_mm_per_rad"],
                    }
                )
                keep_mask.append(False)
            else:
                keep_mask.append(True)
    else:
        keep_mask = [True] * len(paths)
    kept = [path for path, keep in zip(paths, keep_mask) if keep]
    for index, path in enumerate(kept):
        path["id"] = f"{channel['id']}_{index}"
        path["parent_channel"] = channel["id"]
    conflict = False
    if target is not None and not kept and paths:
        # Every fibre pulls the wrong way.  Keep the strongest one and say so,
        # because a silent empty channel would look like a routing problem.
        conflict = True
        best = max(paths, key=lambda item: abs(item["moment_arm_mm_per_rad"]))
        best["id"] = f"{channel['id']}_0"
        best["parent_channel"] = channel["id"]
        kept = [best]
        dropped = [item for item in dropped if item["id"] != best["id"]]

    return {
        "id": channel["id"],
        "group": channel["group"],
        "enabled": bool(kept),
        "elements": list(channel["elements"]),
        "expected_origin_bone": preferred[0],
        "expected_insertion_bone": preferred[1],
        "target_joint": channel["target_joint"],
        "target": (
            None
            if target is None
            else {"joint": target[0], "axis": target[1], "action_sign": target[2]}
        ),
        "channels": len(kept),
        "max_force_n": float(channel["max_force_n"]),
        "travel_m": float(channel["travel_m"]),
        "lattice_muscle": B.LATTICE_MUSCLE,
        "lattice_tendon": B.LATTICE_TENDON,
        "series_stiffness_n_per_m": _series_stiffness(kept),
        "sign_conflict": conflict,
        "dropped_paths": dropped,
        "paths": kept,
    }


def _unique_faces(triangles: np.ndarray) -> np.ndarray:
    """Face indices for a triangle soup, built from its unique vertices."""

    vertices, inverse = np.unique(
        np.round(triangles.reshape(-1, 3), 9), axis=0, return_inverse=True
    )
    return inverse.reshape(-1, 3)


if __name__ == "__main__":
    main()
