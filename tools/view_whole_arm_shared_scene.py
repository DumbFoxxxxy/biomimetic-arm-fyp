"""View the transferred full arm or replay a checked geometric motion."""

import argparse
import sys
import time
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read

sys.path.insert(0, str(ROOT))
from project_runtime import load_mujoco_model, load_mujoco_xml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motion", help="Name in whole_arm_shared_motion.json")
    parser.add_argument("--scapular-action", choices=("protraction", "retraction_alternative", "upward_rotation"),
                        help="Replay a checked scapular action with eight fixed-end wires")
    parser.add_argument("--shoulder-action", choices=("K2_GH_flexion", "K3_GH_abduction"),
                        help="Replay a checked shoulder action with eight fixed-end wires")
    parser.add_argument("--shared-elbow", action="store_true", help="Replay K1 on the fixed K7 assembly")
    parser.add_argument("--cap-relief", action="store_true", help="Use the checked sliding-cap relief in K1 or K7")
    parser.add_argument("--first-pass", action="store_true", help="View the assembled first-pass model at neutral")
    parser.add_argument("--passive-elbow", action="store_true", help="Replay the checked elbow trace with extended rails")
    parser.add_argument("--base", action="store_true", help="Show the checked tabletop base")
    parser.add_argument("--girdle-tpu", action="store_true", help="Show the two new lattices at neutral")
    parser.add_argument("--routed-k7", action="store_true", help="Replay the checked eight-wire portion of K7")
    parser.add_argument("--cuff-support", action="store_true", help="Show the integrated FJ1500 fixed support at neutral")
    parser.add_argument("--combined-cuffs", action="store_true", help="Show both tensioned cuff chords in the routed K7 replay")
    parser.add_argument("--latest-tpu", action="store_true", help="Show the eight meshes from the latest K7 TPU audit at neutral")
    parser.add_argument("--write-assembly", action="store_true", help="Save the latest TPU neutral assembly as MJCF instead of opening the viewer")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    folder = DEST.parent / "whole_arm_first_pass"
    check = read(folder / "whole_arm_shared_scene_check.json")
    path = ROOT / check["model"]["file"]
    if provenance(path) != check["model"]:
        raise ValueError("Shared scene changed since the load check")
    if sum(bool(replay) for replay in (args.passive_elbow, args.motion, args.scapular_action,
                                      args.shoulder_action, args.shared_elbow)) > 1:
        parser.error("Select one replay")
    if args.cap_relief and not (args.shared_elbow or (args.motion == "K7_SC1_pos" and args.combined_cuffs)):
        parser.error("--cap-relief requires --shared-elbow or the combined K7 replay")
    routed_action = args.scapular_action or args.shoulder_action
    is_replay = args.motion or args.passive_elbow or routed_action or args.shared_elbow
    if args.first_pass and (is_replay or args.latest_tpu or args.girdle_tpu or args.cuff_support or args.write_assembly or args.base):
        parser.error("--first-pass opens the complete neutral view and already includes the base")
    if args.latest_tpu and (is_replay or args.girdle_tpu):
        parser.error("--latest-tpu shows the complete neutral TPU geometry")
    if args.write_assembly and not args.latest_tpu:
        parser.error("--write-assembly requires --latest-tpu")
    if args.latest_tpu or routed_action or (args.cap_relief and not args.shared_elbow):
        args.base = True
    if args.routed_k7 and args.motion != "K7_SC1_pos":
        parser.error("--routed-k7 requires --motion K7_SC1_pos")
    if args.girdle_tpu and is_replay:
        parser.error("The girdle TPU mesh overlay shows neutral geometry")
    if args.combined_cuffs and not args.routed_k7:
        parser.error("--combined-cuffs requires --routed-k7")
    if args.cuff_support and is_replay:
        parser.error("--cuff-support is a neutral assembly view")
    has_support = args.cuff_support or args.combined_cuffs or args.latest_tpu or routed_action
    if args.base or args.girdle_tpu or args.passive_elbow or has_support:
        xml = ET.parse(path).getroot()
    if has_support:
        support = read(folder / "whole_arm_infraspinatus_support.json")
        if support["neutral_hits_mm3"]:
            raise ValueError("Infraspinatus support did not pass neutral geometric fit")
        source = support["source"]
        if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
            raise ValueError("Shared infraspinatus support mesh changed")
        old_file = support["replaced_adapter"]["file"]
        asset = next(node for node in xml.findall("./asset/mesh") if node.get("file") == old_file)
        mesh_name = asset.get("name")
        asset.set("file", source["file"])
        offset = np.subtract(source["translation_source_mm"], check["neutral_body_origins_mm"]["humerus"]) * .001
        xml.find(f".//geom[@mesh='{mesh_name}']").set("pos", " ".join(map(str, offset)))
    if args.combined_cuffs or args.latest_tpu or routed_action:
        combined = read(folder / "whole_arm_combined_cuff_k7.json")
        if combined["support"] != provenance(folder / "whole_arm_infraspinatus_support.json"):
            raise ValueError("Combined cuff replay uses another support")
        prior = read(ROOT / support["audit"]["file"])
        revision_path = ROOT / prior["humeral_revision"]["file"]
        if provenance(revision_path) != prior["humeral_revision"]:
            raise ValueError("Cuff hardware reference changed")
        revision = read(revision_path)["selected"]["meshes"]
        transfer = read(ROOT / check["transfer"]["file"])
        for key, source in revision.items():
            if not key.startswith("link:") or key.endswith(("bone_groove_ring", "bone_uncrimped_sleeve")):
                continue
            if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
                raise ValueError(f"Changed cuff hardware: {key}")
            original = transfer["added_sources"][key]
            asset = next(node for node in xml.findall("./asset/mesh") if node.get("file") == original["file"])
            asset.set("file", source["file"])
            offset = np.subtract(source["translation_source_mm"], check["neutral_body_origins_mm"]["humerus"]) * .001
            xml.find(f".//geom[@mesh='{asset.get('name')}']").set("pos", " ".join(map(str, offset)))
    if args.latest_tpu:
        audit = read(folder / "whole_arm_k7_teres_end_relief.json")
        if audit["combined"] != provenance(folder / "whole_arm_combined_cuff_k7.json"):
            raise ValueError("The TPU audit is stale")
        old = xml.find(".//body[@name='neutral_TPU_reference']")
        asset = xml.find("asset")
        old_meshes = {geom.get("mesh") for geom in old.findall("geom")}
        for mesh in list(asset):
            if mesh.get("name") in old_meshes:
                asset.remove(mesh)
        old.clear()
        old.set("name", "neutral_TPU_reference")
        source_map = audit["tpu_sources"]
        for key, source in source_map.items():
            if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
                raise ValueError(f"K7 TPU mesh changed: {key}")
            name = key.replace(":", "_")
            ET.SubElement(asset, "mesh", name=name, file=source["file"], scale=".001 .001 .001")
            ET.SubElement(old, "geom", name=name, type="mesh", mesh=name, rgba=".78 .24 .18 1",
                          pos=" ".join(map(str, np.asarray(source["translation_source_mm"]) * .001)),
                          group="3", contype="0", conaffinity="0")
    if args.passive_elbow:
        audit = read(folder / "whole_arm_passive_elbow_audit.json")
        if audit["scene"] != provenance(folder / "whole_arm_shared_scene_check.json"):
            raise ValueError("The elbow replay was checked on a different scene")
        rail_path = ROOT / audit["extended_rails"]["file"]
        if provenance(rail_path) != audit["extended_rails"]:
            raise ValueError("The extended rail record changed")
        rails = read(rail_path)["shafts"]
        assets = {node.get("file"): node for node in xml.findall("./asset/mesh")}
        for key, rail in rails.items():
            original = audit["nominal_scene_rail_sources"][key]
            if provenance(ROOT / rail["file"])["sha256"] != rail["sha256"]:
                raise ValueError(f"Extended rail mesh changed: {key}")
            node = assets[original["file"]]
            name = node.get("name")
            node.set("file", rail["file"])
            offset = np.subtract(rail["translation_source_mm"], check["neutral_body_origins_mm"]["humerus"]) * .001
            xml.find(f".//geom[@mesh='{name}']").set("pos", " ".join(map(str, offset)))
        flexible_path = folder / "whole_arm_elbow_flexible_audit.json"
        if flexible_path.is_file():
            flexible = read(flexible_path)
            if flexible["replay"] != provenance(folder / "whole_arm_passive_elbow_audit.json"):
                raise ValueError("Elbow retention audit is stale")
            transfer = read(ROOT / check["transfer"]["file"])
            for key, source in flexible["revisions"].items():
                if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
                    raise ValueError(f"Elbow retention mesh changed: {key}")
                original = transfer["added_sources"].get(key)
                if original:
                    node = assets[original["file"]]
                    name = node.get("name")
                    geom = xml.find(f".//geom[@mesh='{name}']")
                    owner = transfer["owner_by_source"][key]
                    node.set("file", source["file"])
                else:
                    owner = "humerus" if key.startswith("series:fixed") else "FJ1477_carriage"
                    name = key.replace(":", "_")
                    ET.SubElement(xml.find("asset"), "mesh", name=name, file=source["file"], scale=".001 .001 .001")
                    body = xml.find(f".//body[@name='{owner}']")
                    geom = ET.SubElement(body, "geom", name=name, mesh=name, type="mesh",
                                         rgba=".53 .63 .72 1", group="1", contype="0", conaffinity="0")
                offset = np.subtract(source["translation_source_mm"],
                                     check["neutral_body_origins_mm"][owner]) * .001
                geom.set("pos", " ".join(map(str, offset)))
    if args.base:
        base = read(folder / "tabletop_base_check.json")
        if (base["scene"] != provenance(folder / "whole_arm_shared_scene_check.json")
                or base["motion"] != provenance(folder / "whole_arm_shared_motion.json")):
            raise ValueError("The base was checked on different geometry or motion")
        source = base["source"]
        if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
            raise ValueError("The tabletop base mesh changed since the fit check")
        if args.motion and base["motion_checks"][args.motion]["first_hit"]:
            raise ValueError("The selected motion intersects the tabletop base")
        ET.SubElement(xml.find("asset"), "mesh", name="tabletop_base",
                      file=source["file"], scale=".001 .001 .001")
        centre = np.asarray(source["translation_source_mm"]) * .001
        ET.SubElement(xml.find(".//body[@name='thorax']"), "geom", name="tabletop_base",
                      type="mesh", mesh="tabletop_base", rgba=".35 .42 .5 1",
                      pos=" ".join(map(str, centre)), contype="0", conaffinity="0", group="1")
    if args.girdle_tpu:
        report = read(folder / "girdle_tpu_geometry.json")
        if report["scene"] != provenance(folder / "whole_arm_shared_scene_check.json"):
            raise ValueError("The TPU meshes were generated on another scene")
        for name, row in report["muscles"].items():
            source = row["source"]
            if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
                raise ValueError(f"The TPU mesh changed: {name}")
            ET.SubElement(xml.find("asset"), "mesh", name=name, file=source["file"], scale=".001 .001 .001")
            centre = np.asarray(source["translation_source_mm"]) * .001
            ET.SubElement(xml.find(".//body[@name='thorax']"), "geom", name=name,
                          type="mesh", mesh=name, rgba=".78 .24 .18 1", group="3",
                          pos=" ".join(map(str, centre)), contype="0", conaffinity="0")
    def use_relief(draft, cap):
        replacement = cap["selected_source"]
        if provenance(ROOT / replacement["file"])["sha256"] != replacement["sha256"]:
            raise ValueError("The relief mesh changed since the geometry audit")
        original_meshes = [mesh for mesh in draft.findall("./asset/mesh") if mesh.get("file") == cap["original"]["file"]]
        if len(original_meshes) != 1 or replacement["translation_source_mm"] != cap["original"]["translation_source_mm"]:
            raise ValueError("The relief does not share the original cap frame")
        original_meshes[0].set("file", replacement["file"])

    if args.cap_relief:
        relief = read(folder / "whole_arm_k1_cap_relief.json")
        if not relief.get("cap_relief", {}).get("geometry_accepted"):
            raise ValueError("The cap relief did not pass the protected-interface checks")
        if not args.shared_elbow:
            if (relief["combined"] != provenance(folder / "whole_arm_combined_cuff_k7.json")
                    or len(relief["K7_relieved_cap"]) != len(combined["combined"])
                    or any(row["hits"] for row in relief["K7_relieved_cap"])):
                raise ValueError("The cap relief did not pass the K7 rigid and wire checks")
            use_relief(xml, relief["cap_relief"])
    if args.first_pass:
        from study_s1_local_force_geometry import load_reference

        first_pass = read(folder / "whole_arm_first_pass_assembly.json")
        first_pass_model = ROOT / first_pass["installed"]["file"]
        if provenance(first_pass_model) != first_pass["installed"]:
            raise ValueError("The first-pass neutral assembly has changed")
        model = load_mujoco_model(first_pass_model)
        if model.ngeom != first_pass["geometry_count"] or np.count_nonzero(model.geom_group == 3) != 8:
            raise ValueError("The first-pass neutral assembly has incomplete geometry")
    elif args.shared_elbow:
        from study_s1_local_force_geometry import load_reference

        k1 = relief if args.cap_relief else read(folder / "whole_arm_k1_shared.json")
        k1_assembly = ROOT / k1["assembly"]["file"]
        if provenance(k1_assembly) != k1["assembly"]:
            raise ValueError("K1 replay uses a different K7 assembly")
        if args.cap_relief:
            draft = ET.parse(k1_assembly).getroot()
            use_relief(draft, k1["cap_relief"])
            model = load_mujoco_xml(ET.tostring(draft, encoding="unicode"), base_dir=k1_assembly.parent)
        else:
            model = load_mujoco_model(k1_assembly)
    elif args.base or args.girdle_tpu or args.passive_elbow or has_support:
        model = load_mujoco_xml(ET.tostring(xml, encoding="unicode"), base_dir=path.parent)
    else:
        model = load_mujoco_model(path)
    data = mujoco.MjData(model)
    if args.latest_tpu:
        if len(source_map) != 8 or np.count_nonzero(model.geom_group == 3) != 8:
            raise ValueError("Latest assembly must have exactly eight audited TPU pieces")
        for name in source_map:
            model.geom(name.replace(":", "_"))
        if args.write_assembly:
            target = folder / "whole_arm_latest_neutral.xml"
            ET.ElementTree(xml).write(target, encoding="utf-8", xml_declaration=True)
            reloaded = load_mujoco_model(target)
            if reloaded.ngeom != model.ngeom or reloaded.nq != model.nq:
                raise ValueError("Exported assembly differs after reloading")
            print(f"Exported and reloaded {reloaded.ngeom} geometries including eight audited TPU meshes")
            return
    timeline = [data.qpos.copy()]
    wire_ports = []
    routed = []
    cuff_frames = []

    def load_wire_overlay(report):
        for item in report["routes"].values():
            points = np.asarray(item["points_source_mm"])
            anchors = [(owner, (points - check["neutral_body_origins_mm"][owner]) * .001)
                       for owner in item["owners"]]
            routed.append((anchors, np.asarray(item["distal_weights"])[:, None]))

    if args.first_pass:
        addresses = load_reference(first_pass["guide_points"])
        load_wire_overlay(load_reference(addresses["routes"]))
    if args.shared_elbow:
        outward = [np.asarray(row["qpos"]) for row in k1["accepted"]]
        if len(outward) < 2 or k1["focal_peak"]["elbow_deg"] <= 0.:
            raise ValueError("K1 has no nonzero checked interval in the K7 assembly")
        timeline = outward + outward[-2::-1]
        k1_combined = load_reference(k1["combined"])
        load_wire_overlay(load_reference(k1_combined["routed"]))

    if routed_action:
        report_name, path_key = (("whole_arm_anatomical_actions.json", "scapular") if args.scapular_action else
                                 ("whole_arm_shoulder_elevation.json", "paths"))
        report = read(folder / report_name)
        if (report["combined"] != provenance(folder / "whole_arm_combined_cuff_k7.json")
                or report["base"] != provenance(folder / "tabletop_base_check.json")):
            raise ValueError("Routed replay was checked against another assembled model")
        path_report = report[path_key].get(routed_action)
        if not path_report or len(path_report["accepted"]) < 2:
            raise ValueError("Selected action has no checked nonzero interval")
        outward = [np.asarray(row["qpos"]) for row in path_report["accepted"]]
        timeline = outward + outward[-2::-1]
        route_path = ROOT / combined["routed"]["file"]
        if provenance(route_path) != combined["routed"]:
            raise ValueError("Routed replay uses another wire geometry")
        load_wire_overlay(read(route_path))
    if args.passive_elbow:
        accepted = audit["samples"][:audit["contiguous_samples"]]
        if not accepted or audit["replacement_neutral_overlaps"]:
            raise ValueError("The elbow variant is blocked at neutral")
        timeline = [np.asarray(row["qpos"]) for row in accepted]
        if audit["first_limit"]:
            timeline += timeline[-2::-1]
    if args.motion:
        report = read(folder / "whole_arm_shared_motion.json")
        if report["model"] != check["model"]:
            raise ValueError("Motion was sampled against a different shared scene")
        if args.motion not in report["paths"]:
            parser.error("Available motions: " + ", ".join(report["paths"]))
        outward = [np.asarray(row["qpos"]) for row in report["paths"][args.motion]["accepted"]]
        if args.routed_k7:
            route_report = read(folder / "whole_arm_posterior_detour_audit.json")
            if route_report["motion"] != provenance(folder / "whole_arm_shared_motion.json"):
                raise ValueError("Routed K7 is stale")
            count = min(row["initial_clear_samples"]
                        for row in route_report["paths"][args.motion]["summary"].values())
            if count < 2:
                raise ValueError("Routed K7 has no nonzero checked samples")
            if args.combined_cuffs:
                if combined["routed"] != provenance(folder / "whole_arm_posterior_detour_audit.json"):
                    raise ValueError("Combined cuff replay uses different routes")
                for record, pose in zip(combined["combined"], report["paths"][args.motion]["accepted"][:count],
                                        strict=True):
                    if record["input_deg"] != pose["input_deg"]:
                        raise ValueError("Combined cuff replay is out of phase")
                    if (record["baseline_wire_hits_mm3"] or record["baseline_wire_pair_hits"]
                            or any(row.get("missing_cuffs") or row["local_cuff_hits"]
                                   or row["cuff_cross_hits_mm3"] or row["cuff_to_neighbor_wires_mm3"]
                                   for row in record["steps"])):
                        break
                    high = next(row for row in record["steps"] if row["load_N"] == 10.)
                    cuff_frames.append([np.asarray(segment) * .001 for segment in
                                        high["central_tendon_tips_world_mm"].values()])
                count = len(cuff_frames)
                if count < 2:
                    raise ValueError("Combined cuff K7 has no nonzero checked samples")
                cuff_frames += cuff_frames[-2::-1]
            outward = outward[:count]
            load_wire_overlay(route_report)
        timeline = outward + outward[-2::-1]
        if not args.routed_k7 and "AC_per_SC" in report["paths"][args.motion]:
            seats = read(folder / "girdle_seat_pair.json")["seats"]
            for route in ("FJ1459_g3", "FJ1536_g2"):
                ports = []
                for role in ("origin", "insertion"):
                    seat = seats[f"{route}_{role}"]
                    owner = seat["owner"]
                    local = np.subtract(seat["port_source_mm"], check["neutral_body_origins_mm"][owner]) * .001
                    ports.append((owner, local))
                wire_ports.append(ports)

    def frame(index):
        data.qpos[:] = timeline[index]
        mujoco.mj_forward(model, data)
        if not np.isfinite(data.geom_xpos).all():
            raise ValueError(f"Nonfinite geometry at frame {index}")

    def draw_wires(scene, frame_index):
        segments = []
        for index, ports in enumerate(wire_ports):
            endpoints = []
            for owner, local in ports:
                body = data.body(owner)
                endpoints.append(body.xmat.reshape(3, 3) @ local + body.xpos)
            color = ([1., .47, .1, 1.] if index == 0 else [.18, .78, 1., 1.])
            segments.append((*endpoints, color))
        for index, (anchors, weights) in enumerate(routed):
            points = []
            for owner, local in anchors:
                body = data.body(owner)
                points.append(local @ body.xmat.reshape(3, 3).T + body.xpos)
            path = (1. - weights) * points[0] + weights * points[1]
            color = [1., .47, .1, 1.] if index % 2 == 0 else [.18, .78, 1., 1.]
            segments.extend((a, b, color) for a, b in zip(path[:-1], path[1:]))
        if cuff_frames:
            segments.extend((a, b, [.63, 1., .30, 1.]) for a, b in cuff_frames[frame_index])
        if len(segments) > scene.maxgeom:
            raise ValueError("Wire overlay exceeds the viewer scene capacity")
        for index, (a, b, color) in enumerate(segments):
            cable = scene.geoms[index]
            mujoco.mjv_initGeom(cable, mujoco.mjtGeom.mjGEOM_CAPSULE, np.zeros(3),
                                np.zeros(3), np.eye(3).ravel(), np.array(color, dtype=np.float32))
            mujoco.mjv_connector(cable, mujoco.mjtGeom.mjGEOM_CAPSULE, .0003, a, b)
        scene.ngeom = len(segments)

    if args.check:
        capacity = len(wire_ports) + sum(len(weights) - 1 for _, weights in routed) + 2 * bool(cuff_frames)
        overlay = mujoco.MjvScene(model, maxgeom=max(1, capacity))
        for index in range(len(timeline)):
            frame(index)
            draw_wires(overlay, index)
            if any(not np.isfinite(overlay.geoms[i].pos).all() for i in range(overlay.ngeom)):
                raise ValueError(f"Nonfinite wire overlay at frame {index}")
        print(f"{model.ngeom} geometries and {len(timeline)} prescribed frames verified")
        return

    from mujoco import viewer
    frame(0)
    with viewer.launch_passive(model, data) as view:
        view.cam.lookat[:] = [-.09, -.08, 1.2]
        view.cam.distance = .84
        view.cam.azimuth = 125
        view.cam.elevation = 8
        # The imported TPU meshes represent neutral geometry only.
        view.opt.geomgroup[3] = int(not is_replay)
        start = time.perf_counter()
        while view.is_running():
            index = int((time.perf_counter() - start) * len(timeline) / 6.) % len(timeline)
            with view.lock():
                frame(index)
                draw_wires(view.user_scn, index)
            view.sync()
            time.sleep(.02)


if __name__ == "__main__":
    main()
