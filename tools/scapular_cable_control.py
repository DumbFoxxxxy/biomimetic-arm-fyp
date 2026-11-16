"""Screen OpenSim-derived local wires and replay cable-driven scapular motion."""

import argparse
import csv
from itertools import combinations, product
import time
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from scipy.spatial import cKDTree

from audit_s1_series_elbow_integration import find_hits, transforms
from audit_whole_arm_k1_shared import closest_gap, section
from build_channel_first_inventory import provenance, read
from build_s1_teres_origin_mount import A, load_component, mesh_geom, rod, save
from compare_opensim_whole_arm import FOLDER, REPORT, ROOT, length, world_points
from project_runtime import load_mujoco_model, load_mujoco_xml
from study_s1_local_force_geometry import load_reference


OUTPUT = FOLDER / "scapular_control"
ASM = FOLDER / "whole_arm_latest_neutral.xml"
RADIUS_MM = .30
RETURN_STIFFNESS = 3.5
RETURN_DAMPING = .52
POINTS = {
    "serratus_middle": "SerratusAnterior_M",
    "serratus_inferior": "SerratusAnterior_I",
    "trapezius_superior": "TrapeziusScapula_S",
    "trapezius_middle": "TrapeziusScapula_M",
    "trapezius_inferior": "TrapeziusScapula_I",
}
AXES = {0: np.array([1., 0., 0.]), 1: np.array([0., 1., 0.]), 2: np.array([0., 0., 1.])}
GAP_CAP_MM = .5


def posed_points(data, origins, route):
    return np.array([world_points(data, item["owner"], np.array(item["port_mm"]) - origins[item["owner"]])
                     for item in route["sites"]])


def make_wire(points):
    return rod(*points, RADIUS_MM) if len(points) == 2 else A.D.outer_envelope(points, RADIUS_MM)


class Geometry:
    def __init__(self):
        combined = read(FOLDER / "whole_arm_combined_cuff_k7.json")
        scene = load_reference(load_reference(combined["motion"])["scene_check"])
        transfer = load_reference(scene["transfer"])
        support = load_reference(combined["support"])
        revision = load_reference(load_reference(support["audit"])["humeral_revision"])["selected"]["meshes"]
        self.sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"],
                        "elbow:FJ1477_bone_adapter": support["source"],
                        "base:T12": read(FOLDER / "tabletop_base_check.json")["source"]}
        self.sources.update({key: entry for key, entry in revision.items()
                             if key.startswith("link:") and key in self.sources})
        self.parts = {key: load_component(entry) for key, entry in self.sources.items()}
        self.owners = {**transfer["owner_by_source"], "base:T12": "thorax"}
        self.origins = {key: np.asarray(value) for key, value in scene["neutral_body_origins_mm"].items()}
        self.routed = load_reference(combined["routed"])["routes"]
        self.model = load_mujoco_model(ASM)
        self.data = mujoco.MjData(self.model)
        self.cached = {}
        self.surfaces = {}
        for owner in ("thorax", "scapula"):
            keys = [key for key in self.parts if key == "bone:scapula" and owner == "scapula" or
                    key.startswith("thorax_reference:") and owner == "thorax"]
            points, ids = [], []
            for key in keys:
                entry = self.sources[key]
                vertices = mesh_geom.load_obj(ROOT / entry["file"])[0] + entry["translation_source_mm"]
                points.append(vertices)
                ids.extend([key] * len(vertices))
            self.surfaces[owner] = (cKDTree(np.concatenate(points)), ids)

    def mount(self, owner, port):
        tree, ids = self.surfaces[owner]
        distance, index = tree.query(port)
        base = tree.data[index]
        return dict(owner=owner, port_mm=np.asarray(port).tolist(), base_mm=base.tolist(),
                    support=ids[index], standoff_mm=float(distance))

    def prescribed(self, action, angle):
        delta = np.zeros(self.model.nv)
        axis, sign, ratio = action["SC_axis"], action["SC_direction"], action["AC_per_SC_same_axis"]
        delta[self.model.joint("SC").dofadr.item() + axis] = sign
        delta[self.model.joint("AC").dofadr.item() + axis] = sign * ratio
        self.data.qpos[:] = self.model.qpos0
        mujoco.mj_integratePos(self.model, self.data.qpos, delta, np.deg2rad(angle))
        mujoco.mj_forward(self.model, self.data)
        return self.data

    def posed_solids(self, data, model=None, excluded=()):
        cache_key = (tuple(data.qpos), tuple(sorted(excluded)))
        if model is None and cache_key in self.cached:
            return self.cached[cache_key]
        model = self.model if model is None else model
        frame = transforms(model, data, self.origins)
        rigid = {key: solid.transform(frame[self.owners[key]]) for key, solid in self.parts.items()}
        wires = {key: section(route, frame)[0] for key, route in self.routed.items() if key not in excluded}
        if model is self.model:
            self.cached[cache_key] = rigid, wires
        return rigid, wires

    def stock(self, route):
        posts = {}
        for index, item in enumerate(route["sites"]):
            # The standoff is a keep-out cylinder with a conforming footprint.
            # The hole and fasteners are left to detailed CAD.
            base, port = map(np.asarray, (item["base_mm"], item["port_mm"]))
            if np.linalg.norm(port - base) < .9:
                return None
            direction = (port - base) / np.linalg.norm(port - base)
            blank = rod(base - .9 * direction, port, .9)
            reference = self.parts[item["support"]]
            if A.I.intersection_volume(blank, reference) <= .001:
                return None
            material = blank - reference
            if len([part for part in material.decompose() if part.volume() > 1e-6]) != 1:
                return None
            posts[f"{route['id']}:post{index}"] = (item["owner"], material)
        return posts

    def screen(self, data, routes, model=None, measure_gap=False):
        candidates = [route for route in routes if not route.get("installed")]
        replaced = {"FJ1459_g3"} if any(r["id"].startswith("serratus") for r in routes) else set()
        rigid, original = self.posed_solids(data, model, replaced)
        frame = transforms(self.model if model is None else model, data, self.origins)
        posed = {r["id"]: posed_points(data, self.origins, r) for r in routes}
        stock = {key: part.transform(frame[owner]) for route in candidates
                 for key, (owner, part) in route["stock"].items()}
        wires = {key: make_wire(points) for key, points in posed.items() if key not in original}
        names = [*rigid, *original, *stock, *wires]
        parts = [*rigid.values(), *original.values(), *stock.values(), *wires.values()]
        old_count = len(rigid) + len(original)
        pairs = np.asarray([(i, j) for i, j in combinations(range(len(names)), 2) if j >= old_count and
                            not _assigned_contact(names[i], names[j], candidates)], dtype=int).reshape(-1, 2)
        hits = find_hits(parts, names, pairs)
        gap = None
        pair = None
        if measure_gap and not hits and len(pairs):
            gap, nearest = closest_gap(parts, pairs, cap=GAP_CAP_MM)
            pair = None if nearest is None else [names[i] for i in nearest]
        return dict(hits=hits, min_nonmating_gap_mm=gap, closest_pair=pair)


def _assigned_contact(left, right, routes):
    # The support blank terminates at the cable eye centre; contact is planned
    # only between that route and its own blanks and at the mated footprints.
    for route in routes:
        for index, item in enumerate(route["sites"]):
            post = f"{route['id']}:post{index}"
            if left == post and right == route["id"]:
                return True
            if left == item["support"] and right == post:
                return True
    return False


def mapped_guides():
    neutral = read(OUTPUT / "opensim_scapula.json")["outbound"][0]["muscles"]
    registrations = read(REPORT)["registration"]
    mapped = {}
    for key, muscle in POINTS.items():
        points = []
        for point in neutral[muscle]["points"]:
            owner = point["frame"]
            mapping = registrations[f"shoulder:{owner}"]
            xyz = np.array(mapping["linear"]) @ point["body_local_mm"] + mapping["translation_mm"]
            points.append((owner, xyz))
        mapped[key] = points[::-1] if key.startswith("serratus") else points
    return mapped


def installed_route(geometry, key):
    route = geometry.routed[key]
    ends = route["points_source_mm"]
    return dict(id=key, installed=True, sites=[dict(owner=owner, port_mm=point)
                                              for owner, point in zip(route["owners"], ends, strict=True)])


def port_options(geometry, owner, source, role):
    tree, _ = geometry.surfaces[owner]
    hints = [np.asarray(source)]
    for offset in ((-12, 0, 0), (-8, -8, 4), (0, 10, 0), (0, 0, 10), (0, -8, -4),
                   (8, 0, 0), (0, 0, -8), (0, 17, 0), (-20, 0, 8), (0, 18, 13), (-18, -10, 12)):
        hints.append(np.asarray(source) + offset)
    candidates = []
    seen = set()
    for hint in hints:
        _, index = tree.query(hint)
        anchor = tree.data[index]
        if role == "wrap":
            directions = [np.array([-1., 0., 0.]), np.array([-1., -1., 0.]), np.array([-1., 1., 0.])]
        else:
            delta = hint - anchor
            directions = ([delta] if np.linalg.norm(delta) > .001 else []) + [
                np.array(direction, dtype=float) for direction in ((-1, 0, 0), (0, -1, 0), (0, 1, 0),
                                                                   (-1, 1, 0), (0, 0, 1), (-1, 0, 1))]
        for direction in directions:
            dr = direction / np.linalg.norm(direction)
            for lift in (3., 6., 10.):
                port = anchor + lift * dr
                item = geometry.mount(owner, port)
                p = tuple(np.round(item["port_mm"], 5))
                if p not in seen and .9 <= item["standoff_mm"] <= 12.:
                    item["reference_distance_mm"] = float(np.linalg.norm(port - source))
                    item["role"] = role
                    seen.add(p)
                    candidates.append(item)
    candidates.sort(key=lambda item: (item["reference_distance_mm"], item["standoff_mm"]))
    return candidates


def candidate_pool(geometry, action, key, points, minimum_leverage_mm):
    site_pools = []
    for index, (owner, point) in enumerate(points):
        role = "origin" if index == 0 else "insertion" if index == len(points) - 1 else "wrap"
        options = port_options(geometry, owner, point, role)
        site_pools.append(options[:22] if role == "wrap" else options[:38])
    angles = [0., 1.5, 3.]
    positions = []
    for angle in angles:
        data = geometry.prescribed(action, angle)
        positions.append([np.array([world_points(data, item["owner"], np.array(item["port_mm"]) -
                                                  geometry.origins[item["owner"]]) for item in options])
                          for options in site_pools])
    indices = np.array(list(product(*(range(len(sites)) for sites in site_pools))))
    lengths = np.zeros((len(angles), len(indices)))
    for index in range(len(points) - 1):
        for row, frame in enumerate(positions):
            lengths[row] += np.linalg.norm(frame[index][indices[:, index]] -
                                           frame[index + 1][indices[:, index + 1]], axis=1)
    arm = -np.diff(lengths, axis=0) / np.deg2rad(1.5)
    selected = np.flatnonzero(arm.min(axis=0) >= minimum_leverage_mm)
    costs = np.sum([[site_pools[i][index]["reference_distance_mm"] for index in indices[:, i]]
                    for i in range(len(points))], axis=0)
    selected = selected[np.argsort(costs[selected], kind="stable")]
    for idx in selected:
        sites = [site_pools[i][j] for i, j in enumerate(indices[idx])]
        yield dict(id=key, sites=sites, installed=False, sampled_min_moment_arm_mm=float(arm[:, idx].min()),
                   summed_reference_distance_mm=float(costs[idx]))


def choose_route(geometry, action, name, points, already=()):
    # A single bearing feature approximates one of the moving anatomy wrap
    # points. Test the two seed choices separately before using more hardware.
    variants = [points] if len(points) == 2 else [[points[0], p, points[-1]] for p in points[1:-1]]
    failures = {}
    for variant in variants:
        for index, route in enumerate(candidate_pool(geometry, action, name, variant, 4.)):
            if index >= 350:
                break
            route["stock"] = geometry.stock(route)
            if route["stock"] is None:
                failures["mount-footprint"] = failures.get("mount-footprint", 0) + 1
                continue
            for angle in (0., 1.5, 3.):
                hits = geometry.screen(geometry.prescribed(action, angle), [*already, route])["hits"]
                if hits:
                    cause = f"{hits[0]['left']} × {hits[0]['right']}"
                    failures[cause] = failures.get(cause, 0) + 1
                    break
            else:
                return route, failures
    return None, failures


def drop_or_hinge(root, name, axis, coefficient, stiffness=None):
    joint = root.find(f".//joint[@name='{name}']")
    parent = next(body for body in root.findall(".//body") if joint in list(body))
    if coefficient == 0.:
        parent.remove(joint)
    else:
        joint.set("type", "hinge")
        joint.set("axis", " ".join(map(str, axis * coefficient)))
        if stiffness is not None:
            joint.set("stiffness", str(stiffness))
            joint.set("damping", str(RETURN_DAMPING))


def active_model(geometry, action, routes):
    root = ET.parse(ASM).getroot()
    root.set("model", "scapular_cable_control")
    axis, sign, ratio = action["SC_axis"], action["SC_direction"], action["AC_per_SC_same_axis"]
    drop_or_hinge(root, "SC", AXES[axis], sign, RETURN_STIFFNESS)
    drop_or_hinge(root, "AC", AXES[axis], float(np.sign(sign * ratio)))
    if ratio:
        equality = ET.SubElement(root, "equality")
        ET.SubElement(equality, "joint", joint1="AC", joint2="SC",
                      polycoef=f"0 {abs(ratio)} 0 0 0", solref=".003 1", solimp=".999 .9999 .0001")
    for frozen in ("GH", "elbow", "FJ1486_swivel", "FJ1477_swivel", "FJ1477_series_travel"):
        drop_or_hinge(root, frozen, AXES[0], 0.)
    tendon, actuator = ET.SubElement(root, "tendon"), ET.SubElement(root, "actuator")
    for route in routes:
        name = route["id"]
        span = ET.SubElement(tendon, "spatial", name=name, width=".0003", rgba=".97 .37 .12 1")
        for index, item in enumerate(route["sites"]):
            owner = item["owner"]
            body = root.find(f".//body[@name='{owner}']")
            pos = (np.array(item["port_mm"]) - geometry.origins[owner]) * .001
            site = f"{name}_{index}"
            ET.SubElement(body, "site", name=site, pos=" ".join(map(str, pos)), size=".001",
                          rgba=".85 .40 .1 1", group="2")
            ET.SubElement(span, "site", site=site)
            if "base_mm" in item:
                base = (np.array(item["base_mm"]) - geometry.origins[owner]) * .001
                ET.SubElement(body, "geom", name=f"{site}_support_envelope", type="capsule", size=".0009",
                              fromto=" ".join(map(str, [*base, *pos])), rgba=".19 .63 .65 .9",
                              group="2", density="0", contype="0", conaffinity="0")
        ET.SubElement(actuator, "motor", name=f"pull_{name}", tendon=name, gear="-1", ctrllimited="true",
                      ctrlrange="0 30", forcelimited="true", forcerange="0 30")
    return root, load_mujoco_xml(ET.tostring(root, encoding="unicode"), base_dir=ASM.parent)


def pose_metrics(data, geometry):
    gh = data.body("humerus").xpos * 1000
    ac = data.body("scapula").xpos * 1000
    line = gh - ac
    reference = geometry.origins["humerus"] - geometry.origins["scapula"]
    angle = np.arctan2(line[2], -line[0]) - np.arctan2(reference[2], -reference[0])
    return dict(coronal_line_deg=float(np.rad2deg(angle)), anterior_mm=float(geometry.origins["humerus"][1] - gh[1]),
                superior_mm=float(gh[2] - geometry.origins["humerus"][2]))


def simulate(model, geometry, action, routes, driven, target_deg, *, neutral_preload=False,
             limits_deg=(-6., 6.), return_brake_Nms=0.):
    data = mujoco.MjData(model)
    def run():
        mujoco.mj_forward(model, data)
        jid = model.joint("SC")
        qadr = jid.qposadr.item()
        vadr = jid.dofadr.item()
        zero_length = data.ten_length.copy()
        eps = 1e-5
        data.qpos[qadr] = eps
        if model.nq == 2:
            data.qpos[model.joint("AC").qposadr.item()] = abs(action["AC_per_SC_same_axis"]) * eps
        mujoco.mj_forward(model, data)
        lever = (zero_length - data.ten_length) / eps
        mujoco.mj_resetData(model, data)
        mujoco.mj_forward(model, data)
        selected = np.array([index for index, route in enumerate(routes) if route["id"] in driven], dtype=int)
        if len(selected) and np.any(lever[selected] <= 0.):
            return dict(usable=False, reason="nonpositive neutral leverage", lever_mm=(1000. * lever).tolist(), frames=[])
        dt = model.opt.timestep
        maximum = 1.8
        return_end = 4.
        frames = []
        peak = 0.
        minimum_angle = 0.
        peak_tensions = np.zeros(model.nu)
        gravity_demand = (float(data.qfrc_bias[vadr]) +
                          abs(action["AC_per_SC_same_axis"]) * float(data.qfrc_bias[1]) if model.nv == 2 else
                          float(data.qfrc_bias[vadr]))
        for step in range(round(return_end / dt) + 1):
            elapsed = step * dt
            current = data.qpos[qadr].item()
            data.ctrl[:] = 0.
            data.qfrc_applied[vadr] = gravity_demand if neutral_preload else 0.
            if elapsed < maximum and len(selected):
                reference = np.deg2rad(target_deg) * min(elapsed / .45, 1.)
                torque = (10. * (reference - current) - .48 * data.qvel[vadr] + RETURN_STIFFNESS * current +
                          (0. if neutral_preload else gravity_demand))
                data.ctrl[selected] = np.clip(torque / len(selected) / lever[selected], 0., 30.)
            elif return_brake_Nms and len(selected):
                # Keep tension during descent to brake motion back toward neutral.
                # The command vanishes as return speed vanishes.
                torque = return_brake_Nms * max(0., -data.qvel[vadr])
                data.ctrl[selected] = np.clip(torque / len(selected) / lever[selected], 0., 30.)
            angle = float(np.rad2deg(current))
            peak = max(peak, angle)
            minimum_angle = min(minimum_angle, angle)
            peak_tensions = np.maximum(peak_tensions, data.ctrl)
            limited = not np.isfinite(data.qpos).all() or not limits_deg[0] <= angle <= limits_deg[1]
            if step % 20 == 0 or step == round(return_end / dt) or limited:
                mujoco.mj_forward(model, data)
                frames.append(dict(time_s=elapsed, qpos=data.qpos.tolist(), sc_deg=float(np.rad2deg(current)),
                                   tension_N=data.actuator_force.tolist(), length_mm=(1000. * data.ten_length).tolist(),
                                   **pose_metrics(data, geometry)))
            if limited:
                return dict(usable=False, reason="exited angular guard", stop_deg=angle, stop_time_s=elapsed,
                            minimum_angle_deg=minimum_angle, peak_deg=peak, peak_tensions_N=peak_tensions.tolist(),
                            lever_mm=(1000. * lever).tolist(), neutral_gravity_demand_Nm=gravity_demand, frames=frames)
            mujoco.mj_step(model, data)
        final = frames[-1]["sc_deg"]
        minimum = 2. if action["SC_axis"] == 1 else 1.
        return dict(usable=peak >= minimum and abs(final) <= .5, minimum_demo_deg=minimum,
                    lever_mm=(1000. * lever).tolist(),
                    neutral_gravity_demand_Nm=gravity_demand, peak_deg=peak, final_deg=final,
                    minimum_angle_deg=minimum_angle, peak_tensions_N=peak_tensions.tolist(), frames=frames)
    return run()


def motion_clearance(geometry, action, routes, stop_deg):
    outbound, returning = [], []
    previous = None
    first = None
    angles = np.linspace(0., stop_deg, int(np.ceil(stop_deg / .25)) + 1)
    for angle in angles:
        data = geometry.prescribed(action, angle)
        probe = geometry.screen(data, routes, measure_gap=True)
        if probe["hits"]:
            first = dict(angle_deg=angle, **probe)
            if previous is not None:
                for mid in np.linspace(previous, angle, 4)[1:-1]:
                    probe = geometry.screen(geometry.prescribed(action, mid), routes, measure_gap=True)
                    if probe["hits"]:
                        first = dict(angle_deg=mid, **probe)
                        break
                    outbound.append(dict(angle_deg=mid, **probe))
            break
        previous = angle
        outbound.append(dict(angle_deg=angle, **probe))
    for row in reversed(outbound):
        probe = geometry.screen(geometry.prescribed(action, row["angle_deg"]), routes)
        if probe["hits"]:
            raise RuntimeError("Cable check differs on return")
        returning.append(row["angle_deg"])
    return dict(outbound=outbound, return_angles_deg=returning, first_limit=first)


def route_loads(data, geometry, routes, result, model):
    peak = {}
    for frame in result["frames"]:
        data.qpos[:] = frame["qpos"]
        mujoco.mj_forward(model, data)
        for route, tension in zip(routes, frame["tension_N"], strict=True):
            points = posed_points(data, geometry.origins, route)
            for index, item in enumerate(route["sites"]):
                name = f"{route['id']}:{index}"
                offsets = [points[j] - points[index] for j in (index - 1, index + 1) if 0 <= j < len(points)]
                resultant = sum(offset / np.linalg.norm(offset) for offset in offsets) * tension
                norm = float(np.linalg.norm(resultant))
                if name not in peak or norm > peak[name]["resultant_N"]:
                    peak[name] = dict(owner=item["owner"], resultant_N=norm,
                                      vector_world_N=resultant.tolist(), at_s=frame["time_s"])
    return peak


def serialize(route):
    return {key: value for key, value in route.items() if key != "stock"}


def replay_clearance(geometry, model, routes, result):
    data = mujoco.MjData(model)
    last = None
    stop = None
    smallest = GAP_CAP_MM
    samples = 0
    for index, frame in enumerate(result["frames"]):
        angle = frame["sc_deg"]
        if last is not None and abs(angle - last) < .2 and index < len(result["frames"]) - 1:
            continue
        data.qpos[:] = frame["qpos"]
        mujoco.mj_forward(model, data)
        checked = geometry.screen(data, routes, model, measure_gap=True)
        samples += 1
        if checked["hits"]:
            stop = dict(frame=index, time_s=frame["time_s"], angle_deg=angle, hits=checked["hits"])
            break
        if checked["min_nonmating_gap_mm"] is not None:
            smallest = min(smallest, checked["min_nonmating_gap_mm"])
        last = angle
    return dict(samples=samples, smallest_new_interface_gap_mm=smallest, gap_cap_mm=GAP_CAP_MM, first_limit=stop)


def full_trial(geometry, actions, key, routes, driven, target, *, suffix=""):
    action = actions[key]
    tree, model = active_model(geometry, action, routes)
    motion = simulate(model, geometry, action, routes, driven, target, neutral_preload=True)
    full_gravity = simulate(model, geometry, action, routes, driven, target)
    no_wire = simulate(model, geometry, action, routes, [], target, neutral_preload=True)
    out = dict(action=key, axis=action["SC_axis"], sign=action["SC_direction"],
               AC_per_SC=action["AC_per_SC_same_axis"], previous_screen_extent_deg=action["accepted"][-1]["input_deg"],
               driven=driven, routes=[route["id"] for route in routes],
               target_deg=target, neutral_preload_case=motion,
               zero_tension_control={field: value for field, value in no_wire.items() if field != "frames"},
               gravity_without_preload={field: value for field, value in full_gravity.items() if field != "frames"})
    print(f"{key}{suffix}: dynamic {motion.get('peak_deg')}, return {motion.get('final_deg')} deg; checking geometry",
          flush=True)
    if motion["frames"]:
        maximum = max(float(frame["sc_deg"]) for frame in motion["frames"])
        out["geometric_pass"] = motion_clearance(geometry, action, routes, max(target + .25, maximum))
        out["dynamic_clearance"] = replay_clearance(geometry, model, routes, motion)
        out["peak_local_reaction"] = route_loads(mujoco.MjData(model), geometry, routes, motion, model)
        out["passed_constrained_demo"] = bool(motion["usable"] and out["geometric_pass"]["first_limit"] is None and
                                              out["dynamic_clearance"]["first_limit"] is None and
                                              no_wire.get("peak_deg", 5.) < .1)
    path = OUTPUT / f"{key}{suffix}.xml"
    path.write_bytes(ET.tostring(tree, encoding="utf-8", xml_declaration=True).replace(b"\r\n", b"\n"))
    out["model"] = provenance(path)
    return out


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    geometry = Geometry()
    action_ref = FOLDER / "whole_arm_anatomical_actions.json"
    actions = read(action_ref)["scapular"]
    k6 = actions["upward_rotation"]
    landmarks = mapped_guides()
    serratus, serratus_fails = choose_route(geometry, k6, "serratus_middle", landmarks["serratus_middle"])
    if not serratus:
        raise RuntimeError(f"No clear serratus route: {serratus_fails}")
    trapezius, trap_fails = choose_route(geometry, k6, "trapezius_superior", landmarks["trapezius_superior"],
                                       already=[serratus])
    if not trapezius:
        raise RuntimeError(f"No clear trapezius route: {trap_fails}")
    rhomboid = installed_route(geometry, "FJ1536_g2")
    previous_serratus = installed_route(geometry, "FJ1459_g3")
    selected = [serratus, trapezius, rhomboid]
    result = {"upward_rotation": full_trial(geometry, actions, "upward_rotation", selected,
                                            [serratus["id"], trapezius["id"]], 2.6)}
    neutral_baselines = {"upward_rotation": ["FJ1536_g2"], "protraction": ["FJ1459_g3"],
                         "retraction_alternative": ["FJ1536_g2"]}
    scoped = {}
    for key, driven in neutral_baselines.items():
        if key == "upward_rotation":
            action = actions[key]
            _, model = active_model(geometry, action, selected)
            scoped["K6_single_line"] = {route: {k: v for k, v in simulate(
                model, geometry, action, selected, [route], 2.6, neutral_preload=True).items() if k != "frames"}
                for route in (serratus["id"], trapezius["id"])}
            continue
        action = actions[key]
        leverages = []
        for route in selected:
            left = length(posed_points(geometry.prescribed(action, 0.), geometry.origins, route))
            right = length(posed_points(geometry.prescribed(action, .5), geometry.origins, route))
            leverages.append((left - right) / np.deg2rad(.5))
        available = [route["id"] for route, task_leverage in zip(selected, leverages, strict=True) if task_leverage > 1.]
        scoped[f"{key}_with_K6_routes"] = dict(route_order=[route["id"] for route in selected],
                                               nominal_leverage_mm=leverages, integrated_followup=bool(available))
        if available:
            candidate = full_trial(geometry, actions, key, selected, available, 1.3, suffix="_k6_routes")
            clearance = candidate.get("dynamic_clearance") or {}
            if candidate["neutral_preload_case"]["usable"] and not clearance.get("first_limit"):
                result[key] = candidate
                continue
            scoped[f"{key}_with_K6_routes"]["failed_trial"] = candidate
        result[key] = full_trial(geometry, actions, key, [previous_serratus, rhomboid], driven, 1.3)
        result[key]["separate_baseline_routes"] = True
    report = dict(assembly=provenance(ASM), source=provenance(OUTPUT / "opensim_scapula.json"),
                  prior_geometry=provenance(action_ref), routes={route["id"]: serialize(route)
                  for route in [*selected, previous_serratus]}, search_failures={serratus["id"]: serratus_fails,
                  trapezius["id"]: trap_fails}, analyses=scoped, trials=result,
                  settings=dict(return_stiffness_Nm_per_rad=RETURN_STIFFNESS,
                                return_damping_Nms_per_rad=RETURN_DAMPING,
                                assumed_tension_limit_N=30., gap_cap_mm=GAP_CAP_MM),
                  scope="Reduced constrained SC/AC MuJoCo motions on original geometry; only positive wire motor input. "
                        "Assumed linear return spring and a separate neutral counterbalance torque. Collision checks "
                        "target new rigid posts and 0.60 mm wires; old rigid geometry uses the matching prior screen. "
                        "Posts are routing support envelopes with a planned eye and mounting interface.")
    save(OUTPUT / "scapular_cable_control.json", report)
    with (OUTPUT / "engineering_anchors.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("route", "index", "owner", "role", "port_X_mm", "port_Y_mm", "port_Z_mm",
                         "base_X_mm", "base_Y_mm", "base_Z_mm", "support", "standoff_mm", "peak_K6_reaction_N"))
        for route in selected:
            for index, item in enumerate(route["sites"]):
                load = result["upward_rotation"]["peak_local_reaction"][f"{route['id']}:{index}"]
                writer.writerow((route["id"], index, item["owner"], item.get("role"), *item["port_mm"],
                                 *item.get("base_mm", [None] * 3), item.get("support"),
                                 item.get("standoff_mm"), load["resultant_N"]))
    for key, item in result.items():
        motion = item["neutral_preload_case"]
        print(key, motion.get("peak_deg"), motion.get("final_deg"), motion.get("usable"),
              (item.get("dynamic_clearance") or {}).get("first_limit"), flush=True)


def viewer(action):
    report = read(OUTPUT / "scapular_cable_control.json")
    trial = report["trials"][action]
    path = ROOT / trial["model"]["file"]
    if provenance(path) != trial["model"]:
        raise ValueError(f"Replay model changed: {path}")
    view_trajectory(path, trial["neutral_preload_case"]["frames"])


def view_trajectory(path, frames):
    model = load_mujoco_model(path)
    data = mujoco.MjData(model)
    from mujoco import viewer as gui
    with gui.launch_passive(model, data) as view:
        view.opt.geomgroup[3] = 0
        view.opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = 1
        view.cam.lookat[:] = [-.075, -.065, 1.29]
        view.cam.distance = .66
        view.cam.azimuth = 145
        view.cam.elevation = 10
        start = time.perf_counter()
        while view.is_running():
            index = int((time.perf_counter() - start) * len(frames) / 5) % len(frames)
            with view.lock():
                data.qpos[:] = frames[index]["qpos"]
                mujoco.mj_forward(model, data)
            view.sync()
            time.sleep(.015)


def check_replays():
    report = read(OUTPUT / "scapular_cable_control.json")
    for key, trial in report["trials"].items():
        path = ROOT / trial["model"]["file"]
        if provenance(path) != trial["model"]:
            raise ValueError(f"Replay model changed: {path}")
        model = load_mujoco_model(path)
        data = mujoco.MjData(model)
        frames = trial["neutral_preload_case"]["frames"]
        for frame in frames:
            data.qpos[:] = frame["qpos"]
            mujoco.mj_forward(model, data)
            np.testing.assert_allclose(1000. * data.ten_length, frame["length_mm"], rtol=0., atol=1e-7)
            if min(frame["tension_N"]) < 0.:
                raise ValueError("Cable motor applied a negative tension")
        interval = [min(row["sc_deg"] for row in frames), max(row["sc_deg"] for row in frames)]
        leverage = []
        for angle in np.linspace(*interval, 21):
            lengths = []
            for perturbation in (-1e-5, 1e-5):
                data.qpos[0] = np.deg2rad(angle) + perturbation
                if model.nq == 2:
                    data.qpos[1] = abs(trial["AC_per_SC"]) * data.qpos[0]
                mujoco.mj_forward(model, data)
                lengths.append(data.ten_length.copy())
            leverage.append(1000. * (lengths[0] - lengths[1]) / 2e-5)
        ranges = np.asarray([np.min(leverage, axis=0), np.max(leverage, axis=0)]).T
        for route, (lower, upper) in zip(trial["routes"], ranges, strict=True):
            if route in trial["driven"] and lower <= 0.:
                raise ValueError(f"{key}: {route} changes pulling sign in the active interval")
            print(f"{key}: {route} checked leverage {lower:.2f} to {upper:.2f} mm")
        print(f'{key}: {len(frames)} replay frames; constrained demo passed: {trial["passed_constrained_demo"]}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--view", choices=("upward_rotation", "protraction", "retraction_alternative"))
    action.add_argument("--check", action="store_true", help="Reload and check the recorded replay models")
    args = parser.parse_args()
    if args.view:
        viewer(args.view)
    elif args.check:
        check_replays()
    else:
        main()
