"""Build lower-profile scapular muscle mounts against the current thorax."""

from itertools import combinations
from pathlib import Path

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, M, export, load_component, rod, save, shoulder_pads
from study_elbow_series_package import bounded, overlaps, union


def connected(part):
    return sum(piece.volume() > 1e-8 for piece in part.as_original().simplify(1e-6).decompose()) == 1


def positive_half(origin, normal, low, high):
    basis = shoulder_pads.pad_basis(normal)
    return A.I._manifold().Manifold.cube([100., 100., high - low]).translate(
        [-50., -50., low]).transform(np.c_[basis, origin])


def material_hits(new, fixed):
    world = bounded(fixed)
    result = {name: hit for name, part in new.items() if (hit := overlaps(part, world))}
    for (name, part), (other, target) in combinations(new.items(), 2):
        amount = float(A.I.intersection_volume(part, target))
        if amount > .001:
            result.setdefault(name, {})[other] = amount
    return result


def low_profile_backing(head, normal, centre, holes, cap_thickness, source_parts, fixed):
    mf = A.I._manifold()
    points = [np.asarray(row["axis_point_source_mm"]) for row in holes]
    rear = [row["bone_crossings"][1]["position_along_cap_normal_mm"] for row in holes]
    bone = fixed["bone:scapula"]
    trials = []
    for outer, radius, width in ((h, r, w) for h in (-10.8, -10.5, -10.2, -9.9, -9.6)
                                 for r in (4., 3.5, 3.1) for w in (2.5, 1.5)):
        if any(outer >= plane for plane in rear):
            continue
        feet = [rod(p + outer * normal, p + (plane + .5) * normal, radius)
                for p, plane in zip(points, rear, strict=True)]
        basis = shoulder_pads.pad_basis(normal)
        p, q = (basis.T @ (point - centre) for point in points)
        span = np.linalg.norm(p[:2] - q[:2])
        bridge = mf.Manifold.cube([width, span, 1.2]).translate(
            [(p[0] + q[0] - width) / 2., min(p[1], q[1]), outer]).transform(np.c_[basis, centre])
        holes_drill = union(rod(point + 4. * normal, point - 35. * normal, 1.25) for point in points)
        backing = (union([*feet, bridge]) - bone) - holes_drill
        candidate = {"backing_plate": backing}
        head_seat = cap_thickness + .75
        for i, point in enumerate(points):
            # Nominal root-core and head. The thread interface is identified by
            # its tapped bore instead of representing overlapping solid teeth.
            tip = outer + .15
            shaft = rod(point + tip * normal, point + (head_seat + .15) * normal, 1.19)
            head_solid = rod(point + head_seat * normal, point + (head_seat + 3.) * normal, 2.75)
            candidate[f"bolt_{i}"] = shaft + head_solid
        row = dict(back_outer_mm=outer, bearing_radius_mm=radius, bridge_width_mm=width,
                   nominal_engagement_at_hole_axes_mm=[plane - outer - .15 for plane in rear])
        trials.append(row)
        if not connected(backing):
            row["disconnected"] = True
            continue
        row["hits_mm3"] = material_hits(candidate, fixed)
        if not row["hits_mm3"]:
            return candidate, trials, row
    return None, trials, None


def reversed_rear_bolts(normal, holes, cap_thickness, fixed):
    mf = A.I._manifold()
    basis = shoulder_pads.pad_basis(normal)
    trials, chosen = [], {}
    front_nut_base = cap_thickness + .65
    nut = (mf.Manifold.cylinder(2.4, 5.5 / np.sqrt(3.), circular_segments=6)
           - mf.Manifold.cylinder(3., 1.6, circular_segments=48).translate([0., 0., -.3]))
    for index, row in enumerate(holes):
        point = np.asarray(row["axis_point_source_mm"])
        rear = row["bone_crossings"][1]["position_along_cap_normal_mm"]
        selected = None
        for stock_depth, head_height, radius in ((d, h, r) for d in (1., .8, .6, .4)
                                                for h in (1.2, .9, .7) for r in (3.5, 3.1, 2.8)):
            shoulder = rear - stock_depth
            foot = rod(point + shoulder * normal, point + (rear + 1.4) * normal, radius)
            foot -= fixed["bone:scapula"] + rod(point + 4. * normal, point - 35. * normal, 1.7)
            head = rod(point + (shoulder - head_height) * normal, point + shoulder * normal, 2.6)
            shaft = rod(point + (shoulder - .1) * normal,
                        point + (front_nut_base + 3.) * normal, 1.5)
            parts = {f"rear_bearing_{index}": foot, f"bolt_{index}": head + shaft,
                     f"front_nut_{index}": nut.transform(np.c_[basis, point + front_nut_base * normal])}
            area = float((head.translate(.1 * normal) ^ foot).volume()) / .1
            trial = dict(hole=index, bearing_axial_depth_at_centre_mm=stock_depth, head_height_mm=head_height,
                         seat_radius_mm=radius, approximate_bearing_area_mm2=area)
            trials.append(trial)
            if not connected(foot) or area <= 0.:
                trial["no_connected_bearing"] = True
                continue
            trial["hits_mm3"] = material_hits(parts, {**fixed, **chosen})
            if not trial["hits_mm3"]:
                chosen.update(parts)
                selected = trial
                break
        if selected is None:
            return None, trials, None
    return chosen, trials, {"style": "rear low heads, front nuts and individual contoured rear bearings"}


def localized_reverse_bolt(normal, hole, cap_thickness, original_plate, fixed):
    mf = A.I._manifold()
    basis = shoulder_pads.pad_basis(normal)
    point = np.asarray(hole["axis_point_source_mm"])
    rear = hole["bone_crossings"][1]["position_along_cap_normal_mm"]
    rib_clearance = []
    for key, bone in fixed.items():
        if key.startswith("thorax_reference:"):
            contact = original_plate ^ bone
            if contact.volume() > .001:
                rib_clearance.append(contact.hull().minkowski_sum(mf.Manifold.cube([.4] * 3, True)))
    relieved = original_plate - union(rib_clearance)
    front_nut_base = cap_thickness + .65
    nut = (mf.Manifold.cylinder(2.4, 5.5 / np.sqrt(3.), circular_segments=6)
           - mf.Manifold.cylinder(3., 1.6, circular_segments=48).translate([0., 0., -.3]))
    nut = nut.transform(np.c_[basis, point + front_nut_base * normal])
    trials = []
    for depth in (1.4, 1.2, 1., .8, .6, .4):
        seat = rear - depth
        pocket = rod(point - 25. * normal, point + seat * normal, 3.2)
        backing = relieved - pocket
        if not connected(backing):
            trials.append(dict(depth_mm=depth, disconnected=True))
            continue
        for height in (1.2, .9, .7):
            for radius in (2.6, 2.35, 2.1):
                head = rod(point + (seat - height) * normal, point + seat * normal, radius)
                shaft = rod(point + (seat - .1) * normal,
                            point + (front_nut_base + 3.) * normal, 1.5)
                area = float((head.translate(.1 * normal) ^ backing).volume()) / .1
                row = dict(depth_mm=depth, head_height_mm=height, head_radius_mm=radius,
                           approximate_bearing_area_mm2=area)
                trials.append(row)
                if area <= 0.:
                    row["no_support"] = True
                    continue
                parts = {"backing_plate": backing, "bolt_1": head + shaft, "front_nut_1": nut}
                row["hits_mm3"] = material_hits(parts, fixed)
                if not row["hits_mm3"]:
                    row["style"] = "original bolt 0; reverse low-head bolt 1 on relieved shared backing"
                    return parts, trials, row
    return None, trials, None


def shortened_front(head, normal, centre, holes, cap_thickness, source_parts, fixed, back_outer):
    points = [np.asarray(row["axis_point_source_mm"]) for row in holes]
    old_support = source_parts["front_support"]
    bone = fixed["bone:scapula"]
    trials = []
    pad_height = 4.
    drill = union(rod(point + 10. * normal, point - 50. * normal, 1.7) for point in points)
    for shift in np.arange(-8., -14.01, -.5):
        pad_top = shift - pad_height
        disc = rod(centre + (pad_top - 2.) * normal, centre + pad_top * normal, 11.)
        lower_support = old_support ^ positive_half(centre, normal, -50., pad_top)
        support = (lower_support + disc) - bone - drill
        candidate = {"cap": source_parts["cap"].translate(normal * shift), "front_support": support}
        if not connected(support):
            trials.append(dict(shift_mm=float(shift), disconnected=True))
            continue
        for bolt_length in (18., 19., 20., 21., 22., 23., 24.):
            candidate = dict(candidate)
            head_seat = cap_thickness + shift + .75
            # Back nut volume is modelled from back_outer-3.15 to -0.75.
            if not back_outer - 4. <= head_seat - bolt_length <= back_outer - 3.15:
                continue
            for i, point in enumerate(points):
                pack = M.m3_parts(point, normal, cap_thickness + shift, back_outer, bolt_length)
                candidate[f"bolt_{i}"] = pack["bolt"]
                candidate[f"front_washer_{i}"] = pack["front_washer"]
            row = dict(shift_mm=float(shift), nominal_bolt_length_mm=bolt_length,
                       support_volume_mm3=float(support.volume()))
            trials.append(row)
            row["hits_mm3"] = material_hits(candidate, fixed)
            if not row["hits_mm3"]:
                return candidate, trials, row
    return None, trials, None


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    parent_file = folder / "assembly_shared_anchors.json"
    parent = read(parent_file)
    sources = parent["sources"]
    rigid = {key: load_component(row) for key, row in sources.items()}
    refs = A.DESIGN / "s1_retained_joint_comparison/shaft_collar/station_60"
    holes_file = refs / "origin_fixation_review/fixation_audit.json"
    mount_file = refs / "origin_mounts/fit_revision/assembly_check.json"
    route_file = A.DESIGN / "s1_port_path_revision/routes_compact_distal_candidate.json"
    routes = {row["head"]: row for row in read(route_file)["routes"]}
    mounts = read(mount_file)["proximal_mounts"]
    changed, dropped, reports = {}, set(), []
    for row in read(holes_file)["rows"]:
        head = row["head"]
        port = routes[head]["ports"]["origin"]
        normal = np.asarray(port["normal_source"])
        centre = np.asarray(port["inner_face_source_mm"])
        cap_thickness = float((np.asarray(port["outer_face_source_mm"]) - centre) @ normal)
        cap_key = f"cap:{head}/origin"
        prefix = f"origin:{head}/"
        existing = {key.removeprefix(prefix): part for key, part in rigid.items() if key.startswith(prefix)}
        existing["cap"] = rigid[cap_key]
        if head == "FJ1500":
            names = {"backing_plate", "bolt_0", "bolt_1", "nut_0", "nut_1",
                     "back_washer_0", "back_washer_1"}
            fixed = {key: part for key, part in rigid.items() if key not in {prefix + name for name in names}}
            candidate, trials, chosen = low_profile_backing(head, normal, centre, row["holes"],
                                                            cap_thickness, existing, fixed)
            reverse_trials = []
            if candidate is None:
                candidate, reverse_trials, chosen = reversed_rear_bolts(normal, row["holes"], cap_thickness, fixed)
            if candidate is None:
                names = {"backing_plate", "bolt_1", "nut_1", "back_washer_1"}
                fixed = {key: part for key, part in rigid.items() if key not in {prefix + name for name in names}}
                candidate, local_trials, chosen = localized_reverse_bolt(
                    normal, row["holes"][1], cap_thickness, existing["backing_plate"], fixed)
            else:
                local_trials = []
        else:
            reverse_trials = []
            local_trials = []
            names = {"front_support", "bolt_0", "bolt_1", "front_washer_0", "front_washer_1"}
            fixed = {key: part for key, part in rigid.items()
                     if key not in {prefix + name for name in names} | {cap_key}}
            candidate, trials, chosen = shortened_front(head, normal, centre, row["holes"], cap_thickness,
                                                        existing, fixed, mounts[head]["back_outer_from_cap_inner_mm"])
        report = dict(head=head, trials=trials, reversed_bolt_trials=reverse_trials,
                      localized_reverse_trials=local_trials, selected=chosen)
        reports.append(report)
        if candidate is None:
            print(head, "no candidate; tried", len(trials), flush=True)
            continue
        roundtrip, replacements = {}, {}
        for name, part in candidate.items():
            key = cap_key if name == "cap" else prefix + name
            for tolerance in (1e-6, 1e-5, 1e-4):
                try:
                    mesh, source = export(part.as_original().simplify(tolerance),
                                          folder / f"thorax_relief_{head}_{name}.obj", centre)
                except (AssertionError, ValueError) as error:
                    report["export_error"] = dict(part=name, error=str(error) or type(error).__name__)
                    continue
                report.pop("export_error", None)
                roundtrip[key], replacements[key] = mesh, source
                break
            else:
                break
        if "export_error" in report:
            print(head, "export failed", report["export_error"], flush=True)
            continue
        report["exported_hits_mm3"] = material_hits(roundtrip, fixed)
        report["sources"] = replacements
        if not report["exported_hits_mm3"]:
            for key in {prefix + name for name in names} - replacements.keys():
                dropped.add(key)
                rigid.pop(key)
            rigid.update(roundtrip)
            changed.update(replacements)
        print(head, "selection", chosen, "readback hits", report["exported_hits_mm3"], flush=True)
    save(folder / "scapular_thorax_relief.json", dict(
        generator=provenance(Path(__file__).resolve()), baseline=provenance(parent_file),
        fixation_geometry=provenance(holes_file), previous_mounts=provenance(mount_file),
        route_reference=provenance(route_file), reports=reports, replacements=changed, replaced_nuts_and_washers=sorted(dropped),
        scope="Neutral rigid candidate. Thorax and scapular bone geometry are unchanged. FJ1500 "
        "uses a shaped tapped backing or rear low-head fasteners and anterior nuts. FJ1504 reuses the existing "
        "bore axes and rear assembly, with a shorter cap stack. TPU pad, local wire, motion, "
        "thread strength, loading and assembly access are separate work."))


if __name__ == "__main__":
    main()
