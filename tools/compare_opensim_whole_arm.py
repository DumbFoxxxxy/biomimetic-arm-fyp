"""Map anatomical references onto the existing whole-arm cable assembly."""

import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.spatial import cKDTree

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, mesh_geom, rod, save
from audit_s1_series_elbow_integration import find_hits, transforms
from audit_whole_arm_k1_shared import CAP, closest_gap, section
from project_runtime import load_mujoco_model, load_mujoco_xml
from study_s1_local_force_geometry import load_reference
from summarize_whole_arm_first_pass import CHANNELS


FOLDER = DEST.parent / "whole_arm_first_pass"
REPORT = FOLDER / "opensim_comparison.json"
INPUT = DEST.parent.parent / "输入"
K1 = FOLDER / "whole_arm_k1_cap_relief.json"
ENGINEERING = "whole_arm_first_pass_assembly.json"
ARM26_SOURCE = "https://github.com/opensim-org/opensim-models/tree/d9b05d470b1a481c222372c85b75772faf8f7792/Models/Arm26"
SHOULDER_SOURCE = "https://github.com/ComputationalBiomechanicsLab/rmr-solver/tree/v1.1/OpenSim%20Models/for%20RMR%20solver"
MODELS = {
    "FJ1486": ("arm26", ("BRA",), "肱骨前面远侧 → 尺骨冠突附近"),
    "FJ1477": ("arm26", ("TRIlat", "TRImed"), "肱骨后面 → 尺骨鹰嘴"),
    "FJ1468": ("shoulder", ("DeltoideusClavicle_A",), "锁骨外侧部 → 肱骨三角肌粗隆"),
    "FJ1467": ("shoulder", ("DeltoideusScapula_M",), "肩峰 → 肱骨三角肌粗隆"),
    "FJ1513": ("shoulder", ("DeltoideusScapula_P",), "肩胛冈 → 肱骨三角肌粗隆"),
    "FJ1507": ("shoulder", ("TeresMajor",), "肩胛下角附近 → 肱骨结节间沟内侧唇"),
    "FJ1459_g3": ("shoulder", ("SerratusAnterior_I", "SerratusAnterior_M", "SerratusAnterior_S"), "外侧肋骨 → 肩胛内侧缘肋面"),
    "FJ1536_g2": ("shoulder", ("Rhomboideus_S", "Rhomboideus_I"), "上胸椎棘突区 → 肩胛内侧缘背面"),
    "FJ1500": ("shoulder", ("Infraspinatus_I", "Infraspinatus_S"), "冈下窝 → 肱骨大结节"),
    "FJ1504": ("shoulder", ("Subscapularis_S", "Subscapularis_M", "Subscapularis_I"), "肩胛下窝 → 肱骨小结节"),
}


def vec(node, tag, factor=1.):
    return np.fromstring(node.findtext(tag), sep=" ") * factor


def body_name(node, tag="socket_parent_frame"):
    return node.findtext(tag).rsplit("/", 1)[-1]


def normalize(vector):
    return vector / np.linalg.norm(vector)


def basis(long_axis, cue):
    u = normalize(long_axis)
    v = normalize(cue - (u @ cue) * u)
    return np.column_stack((u, v, np.cross(u, v)))


def axis_fit(source, target, source_cue, target_cue):
    a, b = np.diff(source, axis=0)[0], np.diff(target, axis=0)[0]
    linear = basis(b, target_cue) @ basis(a, source_cue).T * (np.linalg.norm(b) / np.linalg.norm(a))
    return linear, target[0] - linear @ source[0]


def landmark_fit(source, target):
    a, b = source - source.mean(axis=0), target - target.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    u[:, -1] *= np.linalg.det(u @ vt)
    rotation = u @ vt
    linear = rotation * np.sum(b * (a @ rotation.T)) / np.sum(a * a)
    return linear, target.mean(axis=0) - linear @ source.mean(axis=0)


def apply_fit(points, fit):
    return points @ fit[0].T + fit[1]


def vtp_points(path):
    points = ET.parse(path).find(".//Points/DataArray")
    if points.get("format") != "ascii":
        raise ValueError(f"Expected ASCII VTP coordinates in {path}")
    return np.fromstring(points.text, sep=" ").reshape(-1, 3) * 1000.


def end_centre(points, origin, direction):
    station = (points - origin) @ direction
    return points[station >= np.quantile(station, .92)].mean(axis=0)


def fit_anatomy(arm, shoulder, arm_path, shoulder_path, targets, origins, elbow_axis, sc):
    elbow = arm.find("./JointSet/objects/CustomJoint[@name='r_elbow']")
    arm_elbow = vec(elbow.find("./frames/PhysicalOffsetFrame[@name='r_humerus_offset']"), "translation", 1000.)
    arm_axis = vec(elbow.find("./SpatialTransform/TransformAxis[@name='rotation1']"), "axis")
    elbow_target = origins["forearm"]
    gh = origins["humerus"]
    hand_axis = normalize(elbow_target - gh)
    ulna_geom = vtp_points(arm_path.parent / "Geometry/arm_r_ulna.vtp")
    radius_geom = vtp_points(arm_path.parent / "Geometry/arm_r_radius.vtp")
    target_ulna = targets["ulna"][0]
    target_radius = targets["radius"][0]
    distal = np.array([0., -1., 0.])
    arm_wrist = end_centre(ulna_geom, np.zeros(3), distal)
    target_wrist = end_centre(target_ulna, elbow_target, hand_axis)
    radial = end_centre(radius_geom, np.zeros(3), -distal)
    target_radial = end_centre(target_radius, elbow_target, -hand_axis)

    markers = {node.get("name"): vec(node, "location", 1000.)
               for node in shoulder.findall("./MarkerSet/objects/Marker")}
    shoulder_elbow = next(node for node in shoulder.findall("./JointSet/objects/*") if node.get("name") == "elbow")
    shoulder_elbow = vec(shoulder_elbow.find("./frames/PhysicalOffsetFrame"), "translation", 1000.)
    scapula = targets["scapula"][0]
    ai = scapula[scapula[:, 2].argmin()]
    aa_zone = scapula[scapula[:, 2] >= 1330.]
    aa = aa_zone[aa_zone[:, 0].argmin()]
    ts_zone = scapula[(scapula[:, 2] >= 1290.) & (scapula[:, 2] <= 1320.) & (scapula[:, 1] >= -30.)]
    ts = ts_zone[ts_zone[:, 0].argmax()]

    # Midline landmarks are approximate on the assembled thorax surface.
    sternum, c7, t8, xiphoid = (targets[key][0] for key in ("manubrium", "c7", "t8", "xiphoid"))
    ij = sternum[np.argmin(np.linalg.norm(sternum - [0., -145., 1328.], axis=1))]
    c7_tip = end_centre(c7, np.zeros(3), np.array([0., 1., 0.]))
    t8_tip = end_centre(t8, np.zeros(3), np.array([0., 1., 0.]))
    px = xiphoid[xiphoid[:, 2].argmin()]
    specs = [
        ("arm26", "r_humerus", "humerus", ("GH", "elbow"), [np.zeros(3), arm_elbow], [gh, elbow_target],
         (arm_axis, elbow_axis), arm_path.parent / "Geometry/arm_r_humerus.vtp"),
        ("arm26", "r_ulna_radius_hand", "ulna", ("elbow", "distal ulna"), [np.zeros(3), arm_wrist],
         [elbow_target, target_wrist], (radial, target_radial - elbow_target), arm_path.parent / "Geometry/arm_r_ulna.vtp"),
        ("shoulder", "humerus", "humerus", ("GH", "elbow"), [markers["HumHead_Center"], shoulder_elbow],
         [gh, elbow_target], (markers["EpL"] - markers["EpM"], [-1., 0., 0.]), shoulder_path.parent / "TSM_humerus.vtp"),
        ("shoulder", "scapula", "scapula", ("AI", "AA", "TS"), [markers[k] for k in ("ai", "aa", "ts")],
         [ai, aa, ts], None, shoulder_path.parent / "TSM_scapula.vtp"),
        ("shoulder", "clavicle", "clavicle", ("SC", "AC"), [markers[k] for k in ("sc", "acc")],
         [sc, origins["scapula"]], ([0., 1., 0.], [0., 0., 1.]), shoulder_path.parent / "TSM_clavicle.vtp"),
        ("shoulder", "thorax", "thorax", ("IJ", "C7", "T8", "PX"), [markers[k] for k in ("ij", "c7", "t8", "px")],
         [ij, c7_tip, t8_tip, px], None, shoulder_path.parent / "TSM_thorax.vtp"),
    ]
    fits, results = {}, {}
    for label, body, target, names, src, dst, cues, mesh in specs:
        src, dst = np.array(src), np.array(dst)
        fit = axis_fit(src, dst, *map(np.asarray, cues)) if cues else landmark_fit(src, dst)
        points = vtp_points(mesh)
        if label == "shoulder":
            geom = shoulder.find(f"./BodySet/objects/Body[@name='{body}']/attached_geometry/Mesh")
            points *= vec(geom, "scale_factors")
        distances = cKDTree(targets[target][0]).query(apply_fit(points, fit))[0]
        residual = np.linalg.norm(apply_fit(src, fit) - dst, axis=1)
        fits[label, body] = (target, fit)
        result = dict(model=label, source_body=body, target=target, source_mesh=provenance(mesh),
                      landmark_names=names, source_landmarks_mm=src.tolist(), target_landmarks_mm=dst.tolist(),
                      orientation_cues=None if cues is None else np.asarray(cues).tolist(),
                      landmark_residual_mm=residual.tolist(), mesh_vertex_residual_median_p90_mm=np.quantile(distances, [.5, .9]).tolist(),
                      linear=fit[0].tolist(), translation_mm=fit[1].tolist())
        if body == "scapula":
            joint = shoulder.find("./JointSet/objects/CustomJoint[@name='GlenoHumeral']/frames/PhysicalOffsetFrame")
            if joint is not None:
                result["held_out_GH_error_mm"] = float(np.linalg.norm(apply_fit(vec(joint, "translation", 1000.), fit) - gh))
        results[f"{label}:{body}"] = result
    return fits, results


def collect_targets(transfer):
    names = dict(ulna="尺骨", radius="桡骨", scapula="肩胛骨", humerus="肱骨", clavicle="锁骨")
    files = {key: INPUT / "骨骼" / f"{name}.obj" for key, name in names.items()}
    targets = {key: mesh_geom.load_obj(path) for key, path in files.items()}
    refs = {key: [provenance(path)] for key, path in files.items()}
    thorax = {key.rsplit(":", 1)[1]: entry for key, entry in transfer["shared_sources"].items()
              if key.startswith("thorax_reference:")}
    targets["thorax"] = (np.concatenate([mesh_geom.load_obj(ROOT / item["file"])[0] + item["translation_source_mm"]
                                         for item in thorax.values()]), None)
    refs["thorax"] = list(thorax.values())
    for key, ident in {"manubrium": "FJ3290", "c7": "FJ3172", "t8": "FJ3174", "xiphoid": "FJ3153"}.items():
        file = INPUT / "骨骼/支撑参考" / f"{ident}.obj"
        item = thorax.get(ident, {**provenance(file), "translation_source_mm": [0., 0., 0.]})
        targets[key] = (mesh_geom.load_obj(ROOT / item["file"])[0] + item["translation_source_mm"], None)
        refs[key] = [item]
    return targets, refs


def anatomical_path(label, document, name, fits, targets, origins, muscle_vertices):
    muscle = document.find(f"./ForceSet/objects/*[@name='{name}']")
    geometry = muscle.find("GeometryPath")
    nodes = geometry.findall("./PathPointSet/objects/*")
    points = []
    for index, node in enumerate(nodes):
        body = body_name(node)
        target, fit = fits[label, body]
        owner = "forearm" if target == "ulna" else target
        location = node.find("location")
        if location is None:
            points.append(dict(name=node.get("name"), type=node.tag, source_frame=body,
                               owner=owner, status="position-dependent; omitted from static overlay"))
            continue
        local = vec(node, "location", 1000.)
        mapped = apply_fit(local, fit)
        item = dict(name=node.get("name"), type=node.tag, source_frame=body, owner=owner,
                    source_local_mm=local.tolist(), mapped_local_mm=(mapped - origins[owner]).tolist())
        if node.tag == "ConditionalPathPoint":
            item.update(coordinate=node.findtext("socket_coordinate"), active_range=node.findtext("range"))
        if index in (0, len(nodes) - 1):
            # Use the local anatomy mesh to constrain the surface search and expose disagreement.
            vertices = targets[target][0]
            muscle_distance = cKDTree(muscle_vertices).query(vertices)[0]
            mapped_distance = np.linalg.norm(vertices - mapped, axis=1)
            eligible_region = np.ones(len(vertices), dtype=bool)
            if name == "BRA" and target == "humerus":
                # The brachialis footprint is on the anterior humeral shaft.
                shaft = np.abs(vertices[:, 2] - mapped[2]) <= 25.
                mid_y = .5 * (vertices[shaft, 1].min() + vertices[shaft, 1].max())
                eligible_region = shaft & (vertices[:, 1] <= mid_y)
            near_muscle = eligible_region & (muscle_distance <= muscle_distance[eligible_region].min() + 4.)
            eligible = np.flatnonzero(near_muscle)
            nearest = eligible[mapped_distance[eligible].argmin()]
            surface = vertices[nearest]
            item.update(role="origin" if index == 0 else "insertion", surface_local_mm=(surface - origins[owner]).tolist(),
                        projection_distance_mm=float(mapped_distance[nearest]),
                        surface_to_bodyparts_muscle_mm=float(muscle_distance[nearest]),
                        unconstrained_mesh_vertex_distance_mm=float(mapped_distance.min()))
        points.append(item)
    wraps = []
    for path in geometry.findall("./PathWrapSet/objects/PathWrap"):
        name = path.findtext("wrap_object").rsplit("/", 1)[-1]
        matches = [(body, wrap) for body in document.findall("./BodySet/objects/Body")
                   for wrap in body.findall(f"./WrapObjectSet/objects/*[@name='{name}']")]
        wraps.append(dict(object=name, method=path.findtext("method"), point_range=path.findtext("range"), shapes=[
            dict(body=body.get("name"), type=wrap.tag, translation_m=wrap.findtext("translation"),
                 xyz_body_rotation_rad=wrap.findtext("xyz_body_rotation"), quadrant=wrap.findtext("quadrant"),
                 radius_m=wrap.findtext("radius"), length_m=wrap.findtext("length"), dimensions_m=wrap.findtext("dimensions"))
            for body, wrap in matches]))
    reversed_order = muscle.get("name").startswith("SerratusAnterior")
    if reversed_order:
        points.reverse()
        points[0]["role"], points[-1]["role"] = "origin", "insertion"
    return dict(muscle=muscle.get("name"), points=points, wrap_objects=wraps,
                reversed_from_xml_for_anatomical_origin_to_insertion=reversed_order)


def load_k1():
    checked = read(K1)
    reference = checked["assembly"]
    path = ROOT / reference["file"]
    relief = checked["cap_relief"]
    selected = relief["selected_source"]
    if provenance(path) != reference or provenance(ROOT / selected["file"])["sha256"] != selected["sha256"]:
        raise ValueError("K1 assembly or relief mesh changed")
    xml = ET.parse(path).getroot()
    asset = xml.find(f"./asset/mesh[@file='{relief['original']['file']}']")
    asset.set("file", selected["file"])
    return checked, load_mujoco_xml(ET.tostring(xml, encoding="unicode"), base_dir=path.parent)


def pose_cable(data, route):
    weight = np.asarray(route["distal_owner_blend_weights"])[:, None]
    posed = [world_points(data, owner, points) for owner, points in
             zip(route["owners"], route["centreline_local_mm_by_owner"], strict=True)]
    return (1. - weight) * posed[0] + weight * posed[1]


def world_points(data, owner, points):
    body = data.body(owner)
    return np.asarray(points) @ body.xmat.reshape(3, 3).T + 1000. * body.xpos


def length(points):
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def k1_geometry(checked, origins):
    combined = load_reference(checked["combined"])
    scene = load_reference(load_reference(combined["motion"])["scene_check"])
    transfer = load_reference(scene["transfer"])
    support = load_reference(combined["support"])
    hardware = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
               "elbow:FJ1477_bone_adapter": support["source"], "base:T12": load_reference(checked["base"])["source"]}
    sources.update({key: value for key, value in hardware.items() if key.startswith("link:") and key in sources})
    sources[CAP] = checked["cap_relief"]["selected_source"]
    parts = {key: load_component(value) for key, value in sources.items()}
    owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
    routes = load_reference(combined["routed"])["routes"]
    moving = {"forearm", "FJ1486_distal_pivot", "FJ1477_distal_pivot", "FJ1477_carriage"}
    elbow = {"FJ1486", "FJ1477"}
    names = [*parts, *routes]
    count = len(parts)
    pairs = np.array([(i, j) for i, key in enumerate(names) for j in range(i + 1, len(names))
                      if (j >= count or owners[key] != owners[names[j]]) and
                      (key in elbow or names[j] in elbow or owners.get(key) in moving or owners.get(names[j]) in moving)])
    wire_pairs = pairs[pairs[:, 1] >= count]

    def placed(model, data, only_elbow=False):
        frame = transforms(model, data, origins)
        return ([part.transform(frame[owners[key]]) if not only_elbow or owners[key] in moving else part
                 for key, part in parts.items()]
                + [section(route, frame)[0] for route in routes.values()])

    def check(model, data):
        solids = placed(model, data, only_elbow=True)
        hits = find_hits(solids, names, pairs)
        return hits, 0. if hits else float(closest_gap(solids, wire_pairs)[0])

    def trial(model, data, endpoints):
        solids = placed(model, data)
        result = {}
        for (route, bundle), points in endpoints.items():
            candidate_name = f"anatomical chord {route}/{bundle}"
            # Compare a nominal 0.6 mm wire with every rigid solid and other wire.
            selected = [index for index, key in enumerate(names) if key != route]
            result[route, bundle] = find_hits([*solids, rod(*points, .30)], [*names, candidate_name],
                                              np.array([(index, len(names)) for index in selected]))
        return result

    return check, trial


def replay_k1(checked, model, routes, origins):
    data = mujoco.MjData(model)
    accepted = checked["accepted"]
    if checked["checked_return"] != len(accepted):
        raise ValueError("K1 has incomplete return verification")
    check, trial = k1_geometry(checked, origins)
    extra, outbound = [], [accepted[0]]
    for left, right in zip(accepted, accepted[1:]):
        subdivisions = max(1, int(np.ceil(abs(right["elbow_deg"] - left["elbow_deg"]))))
        for fraction in np.arange(1, subdivisions) / subdivisions:
            data.qpos[:] = (1. - fraction) * np.asarray(left["qpos"]) + fraction * np.asarray(right["qpos"])
            mujoco.mj_forward(model, data)
            hits, gap = check(model, data)
            new = dict(qpos=data.qpos.tolist(), elbow_deg=float(np.rad2deg(data.joint("elbow").qpos.item())),
                       hits=hits, moving_wire_interface_gap_capped_3mm=gap)
            extra.append(new)
            if hits:
                raise ValueError(f"K1 newly sampled collision at {new['elbow_deg']}: {hits}")
            outbound.append(new)
        outbound.append(right)
    delta = np.diff([row["elbow_deg"] for row in outbound])
    keys = ("FJ1486", "FJ1477")
    elbow = model.joint("elbow").dofadr.item()
    duration = 1e-5
    values = {key: [] for key in keys}
    isolated = {key: [] for key in keys}
    forward = [(index, False) for index in range(len(outbound))]
    returning = [(index, True) for index in range(len(outbound) - 1, -1, -1)]
    for index, is_return in forward + returning:
        row = outbound[index]
        data.qpos[:] = row["qpos"]
        mujoco.mj_forward(model, data)
        baseline = {key: length(pose_cable(data, routes[key])) for key in keys}
        if "length_mm" in row:
            for key, measured in baseline.items():
                np.testing.assert_allclose(measured, row["length_mm"][key], atol=1e-6, rtol=0.)
        if is_return:
            if "length_mm" not in row and check(model, data) != (row["hits"], row["moving_wire_interface_gap_capped_3mm"]):
                raise ValueError("Interpolated K1 return differs from outward check")
            continue
        offset = np.zeros(model.nv)
        offset[elbow] = 1.
        mujoco.mj_integratePos(model, data.qpos, offset, duration)
        mujoco.mj_forward(model, data)
        for key, measured in baseline.items():
            values[key].append(measured)
            isolated[key].append(-(length(pose_cable(data, routes[key])) - measured) / duration)
    metrics = {}
    bending = delta > .05
    for key in keys:
        gradient = -np.diff(values[key])[bending] / np.deg2rad(delta[bending])
        metrics[key] = dict(net_mm=values[key][-1] - values[key][0], peak_to_peak_mm=float(np.ptp(values[key])),
                            coupled_leverage_mm=[float(gradient.min()), float(gradient.max())],
                            isolated_leverage_mm=[min(isolated[key]), max(isolated[key])])
    return metrics, outbound, extra, trial


def checked_shoulder_actions():
    anatomical = FOLDER / "whole_arm_anatomical_actions.json"
    elevation = FOLDER / "whole_arm_shoulder_elevation.json"
    combined_path = FOLDER / "whole_arm_combined_cuff_k7.json"
    anatomy, shoulder, combined = map(read, (anatomical, elevation, combined_path))
    for data in (anatomy, shoulder):
        if data["combined"] != provenance(combined_path):
            raise ValueError("Shoulder evidence uses another combined assembly")
    motion = load_reference(combined["motion"])
    ranges, playback = {}, {}
    for items, file in ((shoulder["paths"], elevation), (anatomy["scapular"], anatomical)):
        for name, item in items.items():
            samples = item["accepted"]
            playback[name] = [row["qpos"] for row in samples]
            ranges[name] = dict(source=provenance(file), continuous_input_deg=samples[-1]["input_deg"],
                                last_pose={key: value for key, value in samples[-1].items()
                                           if key.endswith("displacement_mm") or key in ("coronal_scapular_deg", "humerus_deg")},
                                first_limit=item["first_limit"], route_metrics=item["line_summary"],
                                outward=len(samples), checked_return=item["checked_return"],
                                max_step_deg=float(max(abs(np.diff([row["input_deg"] for row in samples])), default=0.)))
    for name, item in anatomy["axial"].items():
        cuff = item["individual_cuff_test"]
        samples = motion["paths"][item["source_path"]]["accepted"]
        cuff_inputs = [row["input_deg"] for row in cuff["accepted"]]
        cuff_extent = abs(cuff_inputs[-1]) if cuff_inputs else None
        playback[name] = [row["qpos"] for row in samples if cuff_extent is not None and abs(row["input_deg"]) <= cuff_extent]
        ranges[name] = dict(source=provenance(anatomical), rigid_range_deg=item["rigid_range_deg"],
                            first_rigid_limit=item["first_rigid_limit"], checked_rigid_return=item["checked_rigid_return"],
                            cuff_samples=cuff["accepted"], route=cuff["id"], first_cuff_limit=cuff["first_limit"],
                            cuff_source=cuff["source"], last_passing_cuff_probe_deg=cuff_extent,
                            cuff_probe_max_step_deg=float(max(abs(np.diff(cuff_inputs)), default=0.)))
    accepted = []
    first = None
    for row in combined["combined"]:
        if (row["baseline_wire_hits_mm3"] or row["baseline_wire_pair_hits"]
                or any(step.get("missing_cuffs") or step["local_cuff_hits"] or step["cuff_cross_hits_mm3"]
                       or step["cuff_to_neighbor_wires_mm3"] for step in row["steps"])):
            first = row
            break
        accepted.append(row)
    k7 = motion["paths"]["K7_SC1_pos"]["accepted"][:len(accepted)]
    np.testing.assert_allclose([row["input_deg"] for row in accepted], [row["input_deg"] for row in k7])
    playback["K7"] = [row["qpos"] for row in k7]
    routed = load_reference(combined["routed"])["paths"]["K7_SC1_pos"]
    ranges["K7"] = dict(source=provenance(combined_path), continuous_input_deg=accepted[-1]["input_deg"],
                        humeral_direction_change_deg=accepted[-1]["humerus_deg"], first_limit=first,
                        fixed_end_first_limits_after_prefix={key: item["first_limit"] for key, item in routed["summary"].items()
                                                              if item["first_limit"] is not None},
                        fixed_end_route_metrics={key: dict(net_mm=routed["samples"][len(k7) - 1]["routes"][key]["length_mm"]
                                                                   - routed["samples"][0]["routes"][key]["length_mm"])
                                                 for key in routed["summary"]},
                        max_step_deg=float(max(abs(np.diff([row["input_deg"] for row in accepted])), default=0.)),
                        loaded_chords_at_10N=[next(step["central_tendon_tips_world_mm"]
                                                   for step in row["steps"] if step["load_N"] == 10.)
                                              for row in accepted],
                        scope="Common ten-route geometric samples; the combined cuff solve was performed outbound.")
    return ranges, playback


def direct_chord_trials(trial, model, data, routes, evidence, playback):
    cases = {"FJ1486": "K1", "FJ1477": "K1", "FJ1459_g3": "upward_rotation",
             "FJ1536_g2": "retraction_alternative"}
    for key, action in cases.items():
        bundle_paths = {ref["muscle"]: ref["points"] for ref in routes[key]["anatomical_paths"]}
        neutral_hits, first_limits, passed = {}, {}, {}
        lengths = {bundle: [] for bundle in bundle_paths}
        accepted_lengths = {bundle: [] for bundle in bundle_paths}
        qposes = playback[action]
        angles = ([float(np.rad2deg(row[model.joint("elbow").qposadr.item()])) for row in qposes] if action == "K1"
                  else [row["input_deg"] for row in load_reference(evidence[action]["source"])["scapular"][action]["accepted"]])
        for index, (qpos, angle) in enumerate(zip(qposes, angles, strict=True)):
            data.qpos[:] = qpos
            mujoco.mj_forward(model, data)
            candidates = {}
            for bundle, points in bundle_paths.items():
                ends = np.vstack([world_points(data, point["owner"], point["surface_local_mm"])
                                  for point in (points[0], points[-1])])
                lengths[bundle].append(length(ends))
                if bundle not in first_limits:
                    candidates[key, bundle] = ends
            if candidates:
                hits = trial(model, data, candidates)
                for (_, bundle), collisions in hits.items():
                    if index == 0:
                        neutral_hits[bundle] = collisions
                    if collisions:
                        first_limits[bundle] = dict(input_deg=angle, hits=collisions)
                    else:
                        passed[bundle] = angle
                        accepted_lengths[bundle].append(lengths[bundle][-1])
        routes[key]["direct_attachment_chord_trials"] = {
            bundle: dict(action=action, neutral_hits=neutral_hits[bundle], first_limit=first_limits.get(bundle),
                         continuous_to_deg=passed.get(bundle),
                         nominal_net_mm_over_engineering_path=lengths[bundle][-1] - lengths[bundle][0],
                         net_mm_over_clear_segment=(accepted_lengths[bundle][-1] - accepted_lengths[bundle][0]
                                                    if accepted_lengths[bundle] else None))
            for bundle in bundle_paths
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm26", type=Path, default=ROOT / "OpenSimDownloads/Arm26/arm26.osim")
    parser.add_argument("--shoulder", type=Path, default=ROOT / "OpenSimDownloads/TSM_subject_noWeight.osim")
    args = parser.parse_args()
    paths = {"arm26": args.arm26.resolve(), "shoulder": args.shoulder.resolve()}
    documents = {key: ET.parse(path).getroot().find("Model") for key, path in paths.items()}
    scene = read(FOLDER / "whole_arm_shared_scene_check.json")
    origins = {key: np.asarray(value) for key, value in scene["neutral_body_origins_mm"].items()}
    targets, target_refs = collect_targets(load_reference(scene["transfer"]))
    checked, model = load_k1()
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    fits, registration = fit_anatomy(documents["arm26"], documents["shoulder"], paths["arm26"], paths["shoulder"],
                                      targets, origins, data.xaxis[model.joint("elbow").id], origins["clavicle"])
    guides_ref = provenance(FOLDER / "whole_arm_guide_points.json")
    guides = load_reference(guides_ref)["data"]
    motion, outbound, extra, trial = replay_k1(checked, model, guides, origins)
    shoulder_checks, playback = checked_shoulder_actions()
    playback["K1"] = [row["qpos"] for row in outbound]
    routes = {}
    table = []
    for key, (label, names, zone) in MODELS.items():
        input_mesh = INPUT / "muscle" / f"{key.split('_')[0]}.obj"
        muscle = mesh_geom.load_obj(input_mesh)[0]
        paths_for_route = [anatomical_path(label, documents[label], name, fits, targets, origins, muscle) for name in names]
        guide = guides[key]
        if "owners" in guide:
            ends = np.asarray(guide["centreline_source_mm"])[[0, -1]]
            owners = guide["owners"]
        else:
            ends = np.array([guide["origin_port_source_mm"], guide["loaded_neutral_tendon_tips_world_mm"][1]])
            owners = [guide["origin_owner"], guide["distal_owner"]]
        for comparison in paths_for_route:
            for point, end in zip((comparison["points"][0], comparison["points"][-1]), ends, strict=True):
                if "surface_local_mm" in point:
                    surface = np.array(point["surface_local_mm"]) + origins[point["owner"]]
                    point["distance_to_engineering_port_mm"] = float(np.linalg.norm(surface - end))
        routes[key] = dict(label=CHANNELS[key][0], anatomical_zone=zone, anatomical_paths=paths_for_route,
                           bodyparts_muscle=provenance(input_mesh), engineering_support=CHANNELS[key][1],
                           engineering_owners=owners, engineering_ports_local_mm=(ends - [origins[o] for o in owners]).tolist(),
                           k1=motion.get(key))
        for comparison in paths_for_route:
            first, last = comparison["points"][0], comparison["points"][-1]
            table.append(dict(route=key, muscle=CHANNELS[key][0], model=label, bundle=comparison["muscle"],
                              anatomical_zone=zone, engineering_support=CHANNELS[key][1],
                              reference_origin_owner=first["owner"], reference_insertion_owner=last["owner"],
                              origin_candidate_local_mm=first.get("surface_local_mm"),
                              insertion_candidate_local_mm=last.get("surface_local_mm"),
                              origin_to_engineering_mm=first.get("distance_to_engineering_port_mm"),
                              insertion_to_engineering_mm=last.get("distance_to_engineering_port_mm")))
    direct_chord_trials(trial, model, data, routes, shoulder_checks, playback)
    for row in table:
        key = row["route"]
        if key in motion:
            row.update(engineering_decision="现有导向及转座保留，人体骨面点作附着参考",
                       motion_case="K1", checked_range_deg=outbound[-1]["elbow_deg"],
                       net_length_mm=motion[key]["net_mm"], pulling_leverage_or_torque_range=motion[key]["isolated_leverage_mm"],
                       metric_unit="mm", first_limit=checked["first_limit"]["hits"][0])
            continue
        cases = {
            "FJ1468": "K2_GH_flexion", "FJ1467": "K3_GH_abduction", "FJ1513": "K3_GH_abduction",
            "FJ1507": "K2_GH_flexion", "FJ1459_g3": "upward_rotation", "FJ1536_g2": "retraction_alternative",
            "FJ1500": "K5_external", "FJ1504": "K4_internal",
        }
        action = cases[key]
        evidence = shoulder_checks[action]
        if "route_metrics" in evidence:
            metrics = evidence["route_metrics"][key]
            limit = (evidence["first_limit"] or {}).get("hits", [])
            row.update(net_length_mm=metrics["net_mm"],
                       pulling_leverage_or_torque_range=metrics["signed_leverage_range_mm"], metric_unit="mm",
                       first_limit=limit[0] if limit else None)
        else:
            torques = [sample["torque_Nm"] for sample in evidence["cuff_samples"]]
            limit = evidence["first_cuff_limit"] or evidence["first_rigid_limit"]
            row.update(net_length_mm=None, pulling_leverage_or_torque_range=[min(torques), max(torques)],
                       metric_unit="N·m at 10 N", first_limit="见 JSON 动作记录" if limit else None)
        row.update(engineering_decision="现有路线保留作工程基线，人体骨面候选待实体校核",
                   motion_case=action,
                   checked_range_deg=evidence.get("continuous_input_deg", evidence.get("last_passing_cuff_probe_deg")))
    report = dict(sources={key: provenance(path) for key, path in paths.items()},
                  source_links=dict(arm26=ARM26_SOURCE, shoulder=SHOULDER_SOURCE),
                  target_meshes=target_refs, registration=registration, guide_points=guides_ref, routes=routes,
                  engineering=provenance(FOLDER / ENGINEERING), shoulder_motion=shoulder_checks, playback=playback,
                  K1=dict(source=provenance(K1), accepted_outbound_count=len(outbound),
                          checked_return=len(outbound), added_geometry_samples=extra,
                          max_step_deg=float(np.max(np.abs(np.diff([row["elbow_deg"] for row in outbound])))),
                          continuous_to_deg=outbound[-1]["elbow_deg"],
                          first_limit=checked["first_limit"],
                          minimum_gap_mm=min(row["moving_wire_interface_gap_capped_3mm"] for row in outbound)),
                  arm26_triceps_long_note="TRIlong starts on the scapular part of Arm26 base and crosses the shoulder. The existing humeral carriage is compared to the lateral and medial heads.")
    save(REPORT, report)
    with REPORT.with_suffix(".csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    print(f"Compared {len(routes)} engineering routes with {len(table)} OpenSim bundles")


if __name__ == "__main__":
    main()
