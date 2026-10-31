"""Derive anatomical cable attachments from the BodyParts3D muscle meshes.

BodyParts3D stores every organ as a separate OBJ element.  For a muscle the
mesh is the segmented muscle volume, so the vertices that sit on a bone are the
anatomical attachment footprint of that muscle on that bone.  This tool

1. checks that the BodyParts3D bones line up with the SolidWorks bone exports,
2. measures the attachment footprint of every shoulder/elbow muscle on every
   candidate bone,
3. writes the footprints in the bone's own part coordinates, ready to be used
   as `insertion_part_m` / guide anchors in the cable manifest.

Nothing here modifies the model; it is a measurement pass.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from prepare_cable_drive_exports import load_part_probe, stl_part_offset

BP3D_DIR = Path(r"D:\XJTLU\isa_BP3D_4.0_obj_99")
BP3D_INDEX = BP3D_DIR / "isa_element_parts.txt"
SW_EXPORT_DIR = Path(r"D:\XJTLU\PSP\model\exports")
SW_BONE_DIR = SW_EXPORT_DIR / "bones_local"
PART_PROBE = SW_EXPORT_DIR / "part_probe.txt"

# Right-hand side of the model.  The non-"M" element files are the right side.
BONES = {
    "clavicle": "clavicle",
    "scapula": "scapula",
    "humerus": "humerus",
    "ulna": "ulna",
}

# Explicit right-side names; the bare concept names cover both sides.
RIGHT_BONE_NAMES = {
    "clavicle": r"^right clavicle$",
    "scapula": r"^right scapula$",
    "humerus": r"^right humerus$",
    "ulna": r"^right ulna$",
    "radius": r"^right radius$",
}

# name regex -> which bone the muscle is expected to attach to, and how.
MUSCLES = {
    "supraspinatus": {"origin": "scapula", "insertion": "humerus"},
    "infraspinatus": {"origin": "scapula", "insertion": "humerus"},
    "subscapularis": {"origin": "scapula", "insertion": "humerus"},
    "teres_minor": {"origin": "scapula", "insertion": "humerus"},
    "teres_major": {"origin": "scapula", "insertion": "humerus"},
    "deltoid_clavicular": {"origin": "clavicle", "insertion": "humerus"},
    "deltoid_acromial": {"origin": "scapula", "insertion": "humerus"},
    "deltoid_spinal": {"origin": "scapula", "insertion": "humerus"},
    "coracobrachialis": {"origin": "scapula", "insertion": "humerus"},
    "pectoralis_major_clavicular": {"origin": "clavicle", "insertion": "humerus"},
    "biceps_long_head": {"origin": "scapula", "insertion": "radius"},
    "biceps_short_head": {"origin": "scapula", "insertion": "radius"},
    "brachialis": {"origin": "humerus", "insertion": "ulna"},
    "triceps_long_head": {"origin": "scapula", "insertion": "ulna"},
    "triceps_lateral_head": {"origin": "humerus", "insertion": "ulna"},
    "triceps_medial_head": {"origin": "humerus", "insertion": "ulna"},
}

MUSCLE_PATTERNS = {
    "supraspinatus": r"^supraspinatus$",
    "infraspinatus": r"^infraspinatus$",
    "subscapularis": r"^subscapularis$",
    "teres_minor": r"^teres minor$",
    "teres_major": r"^teres major$",
    "deltoid_clavicular": r"^clavicular part of deltoid$",
    "deltoid_acromial": r"^acromial part of deltoid$",
    "deltoid_spinal": r"^spinal part of deltoid$",
    "coracobrachialis": r"^coracobrachialis$",
    "pectoralis_major_clavicular": r"^clavicular part of pectoralis major$",
    "biceps_long_head": r"^long head of biceps brachii$",
    "biceps_short_head": r"^short head of biceps brachii$",
    "brachialis": r"^brachialis$",
    "triceps_long_head": r"^long head of triceps brachii$",
    "triceps_lateral_head": r"^lateral head of triceps brachii$",
    "triceps_medial_head": r"^medial head of triceps brachii$",
}


def load_index(path: Path) -> list[dict[str, str]]:
    rows = []
    for line in path.read_text(encoding="utf-8-sig").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        rows.append({"fma": parts[0], "name": parts[1].strip(), "file": parts[2].strip()})
    return rows


def element_path(file_id: str) -> Path:
    return BP3D_DIR / f"{file_id}.obj"


def right_side_files(rows: list[dict[str, str]], pattern: str) -> list[str]:
    """Element files whose part name matches and which are not the left side."""
    hits = []
    for row in rows:
        name = row["name"]
        if not re.search(pattern, name):
            continue
        if name.startswith("left ") or name.endswith(" of left") or " left " in name:
            continue
        hits.append(row["file"])
    out = []
    for file_id in hits:
        if element_path(file_id).exists() and file_id not in out:
            out.append(file_id)
    return out


def read_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    vertices: list[list[float]] = []
    faces: list[list[int]] = []
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                vertices.append([float(x), float(y), float(z)])
            elif line.startswith("f "):
                tokens = line.split()[1:]
                idx = [int(token.split("/")[0]) - 1 for token in tokens]
                for k in range(1, len(idx) - 1):
                    faces.append([idx[0], idx[k], idx[k + 1]])
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=np.int64)


def read_stl(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if data[:5].lower() != b"solid" or b"facet" not in data[:2048].lower():
        count = int.from_bytes(data[80:84], "little")
        record = np.dtype([("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
        tri = np.frombuffer(data[84 : 84 + count * 50], dtype=record)["v"].astype(float)
        return tri
    vertices = []
    for line in data.decode("utf-8", errors="ignore").splitlines():
        line = line.strip()
        if line.startswith("vertex"):
            vertices.append([float(v) for v in line.split()[1:4]])
    arr = np.asarray(vertices, dtype=float)
    return arr.reshape(-1, 3, 3)


def sample_surface(tri: np.ndarray, per_triangle: int = 4) -> np.ndarray:
    """Dense-ish surface samples: centroid plus edge midpoints to start with."""
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    samples = [(a + b + c) / 3.0]
    if per_triangle >= 4:
        samples += [(a + b) / 2.0, (b + c) / 2.0, (c + a) / 2.0]
    if per_triangle >= 7:
        samples += [a, b, c]
    return np.vstack(samples)


def load_bp3d_bone(name: str, rows: list[dict[str, str]]) -> np.ndarray:
    files = right_side_files(rows, RIGHT_BONE_NAMES[name])
    if not files:
        raise FileNotFoundError(f"No BodyParts3D element for bone {name}")
    tri = []
    for file_id in files:
        verts, faces = read_obj(element_path(file_id))
        tri.append(verts[faces])
    return np.vstack(tri)


def load_sw_bone(name: str, probe: dict[str, dict]) -> np.ndarray:
    """SolidWorks bone mesh in part coordinates (millimetres).

    The STL exporter writes the mesh in its own output frame, which is offset
    from the part origin.  The offset recovered from `part_probe.txt` is
    exactly what `prepare_cable_drive_exports.py` applies, so adding it here
    puts the attachment points in the same frame as `insertion_part_m`.
    """
    candidates = [SW_BONE_DIR / f"{name}.STL", SW_BONE_DIR / f"{name}.stl"]
    for candidate in candidates:
        if candidate.exists():
            mesh_mm = read_stl(candidate)
            entry = probe.get(name)
            if entry:
                try:
                    offset_mm = stl_part_offset(candidate, entry) * 1000.0
                except ValueError as error:
                    print(f"  ! {name}: offset check failed ({error}); using raw mesh")
                    offset_mm = np.zeros(3)
                mesh_mm = mesh_mm + offset_mm
            return mesh_mm
    raise FileNotFoundError(f"Missing SolidWorks bone export for {name}")


def bbox(points: np.ndarray) -> list[float]:
    return [float(v) for v in np.round(points.max(axis=0) - points.min(axis=0), 2)]


def end_regions(
    vertices: np.ndarray, fraction: float = 0.12
) -> tuple[np.ndarray, np.ndarray]:
    """Split a muscle mesh into the two ends where it meets its bones.

    A muscle is a spindle: its attachments are at the two ends along its
    principal axis.  Proximity alone is not enough, because a broad muscle such
    as the deltoid lies against the humerus all the way up to the shoulder
    joint, and its true insertion is the distal end.
    """

    centred = vertices - vertices.mean(axis=0)
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    projection = centred @ vectors[:, -1]
    low, high = np.percentile(projection, [fraction * 100.0, 100.0 - fraction * 100.0])
    return vertices[projection <= low], vertices[projection >= high]


def attach_ends(
    vertices: np.ndarray,
    fits: dict[str, dict],
    trees: dict[str, cKDTree],
) -> list[dict]:
    """Which bone each end of the muscle sits on, and where."""

    found = []
    for label, region in zip(("end_a", "end_b"), end_regions(vertices)):
        centroid = region.mean(axis=0)
        best: dict | None = None
        for bone, fit in fits.items():
            local = fit["scale"] * (centroid @ fit["rotation"].T) + fit["translation"]
            distance, index = trees[bone].query(local)
            candidate = {
                "label": label,
                "bone": bone,
                "vertices": int(len(region)),
                "centre_part_mm": [round(float(v), 2) for v in local],
                "surface_part_mm": [round(float(v), 2) for v in trees[bone].data[index]],
                "distance_to_bone_mm": round(float(distance), 2),
            }
            if best is None or distance < best["distance_to_bone_mm"]:
                best = candidate
        found.append(best)
    return found


def pca_frame(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centre = points.mean(axis=0)
    centred = points - centre
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    order = np.argsort(values)[::-1]
    return centre, vectors[:, order]


def principal_length(points: np.ndarray) -> float:
    """Length of the point cloud along its principal axis."""
    centred = points - points.mean(axis=0)
    values, _ = np.linalg.eigh(centred.T @ centred / len(centred))
    return float(2.0 * np.sqrt(values[-1]) * np.sqrt(3.0))


def principal_extent(points: np.ndarray) -> float:
    """Peak-to-peak span of the point cloud along its principal axis."""
    centred = points - points.mean(axis=0)
    values, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    axis = vectors[:, -1]
    projection = centred @ axis
    return float(projection.max() - projection.min())


def axis_extents(points: np.ndarray) -> np.ndarray:
    """Peak-to-peak spans along the three principal axes, longest first."""
    centred = points - points.mean(axis=0)
    _, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    axes = vectors[:, ::-1].T
    return np.array(
        [np.ptp(centred @ axis) for axis in axes],
        dtype=float,
    )


def robust_extent(points: np.ndarray, low: float = 2.0, high: float = 98.0) -> float:
    """Principal span between two percentiles, so stray vertices cannot inflate it."""
    centred = points - points.mean(axis=0)
    _, vectors = np.linalg.eigh(centred.T @ centred / len(centred))
    projection = centred @ vectors[:, -1]
    return float(
        np.percentile(projection, high) - np.percentile(projection, low)
    )


def mesh_volume(tri: np.ndarray) -> float:
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    return float(abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return s, R, t so that s * R @ source + t best matches target."""
    source_centre = source.mean(axis=0)
    target_centre = target.mean(axis=0)
    covariance = (source - source_centre).T @ (target - target_centre)
    u, singular, vt = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, sign])
    rotation = vt.T @ correction @ u.T
    scale = float((singular * np.array([1.0, 1.0, sign])).sum())
    scale /= float(((source - source_centre) ** 2).sum())
    translation = target_centre - scale * rotation @ source_centre
    return scale, rotation, translation


def register_bone(
    source: np.ndarray, target_tree: cKDTree, with_scale: bool, iterations: int = 60
) -> dict:
    """Register a BodyParts3D bone onto a SolidWorks bone mesh.

    A thin shell makes the scale badly identifiable: shrinking a flat bone by
    half still leaves its surface within a few millimetres of the target, so a
    free scale search happily collapses.  The scale therefore comes from the
    principal extent ratio and only the pose is optimised.  Quality is scored
    bidirectionally so a shrunken match cannot win.
    """
    source = source.astype(float)
    step = max(1, len(source) // 4000)
    sample = source[::step]
    source_centre, source_frame = pca_frame(sample)
    target_centre, target_frame = pca_frame(
        target_tree.data[:: max(1, len(target_tree.data) // 4000)]
    )
    # The raw principal span is inflated by stray vertices in some exported
    # meshes, so the size prior uses percentile bounds instead.
    scale = robust_extent(target_tree.data) / robust_extent(sample)
    if with_scale:
        tolerance = 0.05
    else:
        tolerance = 0.0

    best: dict | None = None
    for signs in itertools.product([1.0, -1.0], repeat=3):
        if abs(np.linalg.det(np.diag(signs)) - 1.0) > 1e-9:
            continue
        base_rotation = target_frame @ np.diag(signs) @ source_frame.T
        for multiplier in (1.0, 1.0 - tolerance, 1.0 + tolerance):
            scale_now = scale * multiplier
            rotation = base_rotation
            translation = target_centre - scale_now * rotation @ source_centre
            for _ in range(iterations):
                moved = scale_now * (sample @ rotation.T) + translation
                distance, index = target_tree.query(moved, workers=-1)
                _, rotation, _ = kabsch(sample, target_tree.data[index])
                translation = (
                    target_tree.data[index].mean(axis=0) - scale_now * rotation @ sample.mean(axis=0)
                )
            moved_full = scale_now * (source @ rotation.T) + translation
            distance, _ = target_tree.query(moved_full, workers=-1)
            reverse_tree = cKDTree(moved_full[:: max(1, len(moved_full) // 4000)])
            reverse, _ = reverse_tree.query(
                target_tree.data[:: max(1, len(target_tree.data) // 4000)], workers=-1
            )
            score = max(
                float(np.percentile(distance, 90)), float(np.percentile(reverse, 90))
            )
            candidate = {
                "scale": scale_now,
                "rotation": rotation,
                "translation": translation,
                "median_mm": round(float(np.median(distance)), 3),
                "p90_mm": round(float(np.percentile(distance, 90)), 3),
                "reverse_p90_mm": round(float(np.percentile(reverse, 90)), 3),
                "score": score,
            }
            if best is None or score < best["score"]:
                best = candidate
    assert best is not None
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band-mm", type=float, default=4.0)
    parser.add_argument("--with-scale", action="store_true")
    parser.add_argument(
        "--report",
        type=Path,
        default=Path(r"D:\XJTLU\PSP\mujoco_arm\exports\cable_drive\bp3d_attachments.json"),
    )
    args = parser.parse_args()

    rows = load_index(BP3D_INDEX)
    print(f"BodyParts3D index: {len(rows)} rows")

    report: dict[str, object] = {
        "band_mm": args.band_mm,
        "similarity_fit": bool(args.with_scale),
        "bones": {},
        "muscles": {},
    }

    bp3d_bones: dict[str, np.ndarray] = {}
    trees: dict[str, cKDTree] = {}
    fits: dict[str, dict] = {}
    sw_tris: dict[str, np.ndarray] = {}
    probe = load_part_probe(PART_PROBE)
    for bone in BONES:
        bp3d_bones[bone] = load_bp3d_bone(bone, rows)
        sw_tri = load_sw_bone(bone, probe)
        sw_tris[bone] = sw_tri
        trees[bone] = cKDTree(sample_surface(sw_tri, 7))
        fits[bone] = register_bone(
            bp3d_bones[bone].reshape(-1, 3), trees[bone], bool(args.with_scale)
        )

    print("\n== bone registration: BodyParts3D bone fitted onto SolidWorks part ==")
    for bone, fit in fits.items():
        bp3d_points = bp3d_bones[bone].reshape(-1, 3)
        sw_points = sw_tris[bone].reshape(-1, 3)
        length_ratio = principal_length(sw_points) / principal_length(bp3d_points)
        extent_ratio = principal_extent(sw_points) / principal_extent(bp3d_points)
        # Scale-free shape descriptor: the three principal spans relative to the
        # longest one.  Equal descriptors mean the same shape up to a uniform
        # scale; different ones mean the shape itself was redrawn.
        bp3d_axes = axis_extents(bp3d_points)
        sw_axes = axis_extents(sw_points)
        bp3d_ratio = bp3d_axes / bp3d_axes[0]
        sw_ratio = sw_axes / sw_axes[0]
        volume_ratio = (mesh_volume(sw_tris[bone]) / mesh_volume(bp3d_bones[bone])) ** (
            1.0 / 3.0
        )
        stats = {
            "bp3d_vertices": int(bp3d_points.shape[0]),
            "shape_ratio": {
                "principal_length": round(float(length_ratio), 4),
                "principal_extent": round(float(extent_ratio), 4),
                "volume_cube_root": round(float(volume_ratio), 4),
            },
            "axis_ratio": {
                "bp3d": [round(float(v), 3) for v in bp3d_ratio],
                "solidworks": [round(float(v), 3) for v in sw_ratio],
            },
            "robust_extent_mm": {
                "bp3d": round(robust_extent(bp3d_points), 1),
                "solidworks": round(robust_extent(sw_points), 1),
            },
            "fit": {
                "scale": round(float(fit["scale"]), 4),
                "translation_mm": [round(float(v), 2) for v in fit["translation"]],
                "median_mm": fit["median_mm"],
                "p90_mm": fit["p90_mm"],
                "reverse_p90_mm": fit.get("reverse_p90_mm"),
            },
        }
        report["bones"][bone] = stats
        print(
            f"  {bone:9s} fitted scale={fit['scale']:6.4f}  "
            f"extent ratio={extent_ratio:6.4f}  "
            f"bp3d extent={principal_extent(bp3d_points):6.1f} mm  "
            f"sw extent={principal_extent(sw_points):6.1f} mm  "
            f"residual median={fit['median_mm']:6.2f} mm  p90={fit['p90_mm']:6.2f} mm  "
            f"reverse p90={fit.get('reverse_p90_mm', float('nan')):6.2f} mm"
        )
        print(
            f"  {'':9s} principal spans (normalised)  bp3d={np.round(bp3d_ratio, 3).tolist()}  "
            f"solidworks={np.round(sw_ratio, 3).tolist()}"
        )
        print(
            f"  {'':9s} robust span (2-98 pct)  bp3d={robust_extent(bp3d_points):6.1f} mm  "
            f"solidworks={robust_extent(sw_points):6.1f} mm"
        )

    print(f"\n== muscle attachment footprints on the SolidWorks bones (band {args.band_mm} mm) ==")
    for muscle, spec in MUSCLES.items():
        files = right_side_files(rows, MUSCLE_PATTERNS[muscle])
        if not files:
            continue
        verts = np.vstack([read_obj(element_path(f))[0] for f in files])
        entry: dict[str, object] = {"files": files, "vertices": int(verts.shape[0]), "ends": {}}
        for role, bone in spec.items():
            if bone not in trees:
                entry["ends"][role] = {"bone": bone, "ok": False, "reason": "bone not in model"}
                continue
            fit = fits[bone]
            local = fit["scale"] * (verts @ fit["rotation"].T) + fit["translation"]
            dist, _ = trees[bone].query(local, workers=-1)
            mask = dist <= args.band_mm
            count = int(mask.sum())
            if count < 5:
                entry["ends"][role] = {"bone": bone, "vertices_touching": count, "ok": False}
                continue
            footprint = local[mask]
            centre = footprint.mean(axis=0)
            # Nearest SolidWorks bone surface point to the footprint centre.
            _, index = trees[bone].query(centre)
            surface = trees[bone].data[index]
            spread = float(np.linalg.norm(footprint - centre, axis=1).max())
            entry["ends"][role] = {
                "bone": bone,
                "vertices_touching": count,
                "share_of_muscle": round(count / verts.shape[0], 3),
                "centre_part_mm": [round(float(v), 2) for v in centre],
                "centre_part_m": [round(float(v) / 1000.0, 6) for v in centre],
                "snap_to_bone_mm": [round(float(v), 2) for v in surface],
                "footprint_radius_mm": round(spread, 2),
                "distance_to_bone_mm": round(float(np.linalg.norm(centre - surface)), 2),
                "ok": True,
            }
        report["muscles"][muscle] = entry
        origin = entry["ends"].get("origin", {})
        insertion = entry["ends"].get("insertion", {})

        def fmt(end: dict) -> str:
            if not end.get("ok"):
                return f"{end.get('bone', '?'):9s} none"
            return (
                f"{end['bone']:9s} n={end['vertices_touching']:5d} "
                f"centre={end['centre_part_mm']} r={end['footprint_radius_mm']:6.2f} mm"
            )

        print(f"  {muscle:26s} origin  -> {fmt(origin)}")
        print(f"  {'':26s} insert  -> {fmt(insertion)}")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nreport: {args.report}")


if __name__ == "__main__":
    main()
