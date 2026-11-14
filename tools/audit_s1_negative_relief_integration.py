"""Recheck the FJ1504 relief in the full rigid shoulder/elbow candidate."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from study_s1_coupled_scapula import ShoulderAssembly
from study_s1_local_force_geometry import load_reference


def trace(model, offsets, sc_per_degree):
    frames = []
    for angle in offsets:
        sc = angle * sc_per_degree
        hits, limit = model.evaluate(angle, sc)
        frames.append(dict(**model.record(angle, sc), hits=hits, limit=limit))
        if hits or limit:
            break
    return frames


def coordinated_segment(model, start, end, max_step=.05):
    a0, s0 = start
    a1, s1 = end
    count = max(1, int(np.ceil((abs(a1 - a0) + np.linalg.norm(s1 - s0)) / max_step)))
    rows = []
    for fraction in np.linspace(0., 1., count + 1):
        angle = a0 + fraction * (a1 - a0)
        sc = s0 + fraction * (s1 - s0)
        hits, limit = model.evaluate(angle, sc)
        rows.append(dict(**model.record(angle, sc), hits=hits, limit=limit))
        if hits or limit:
            break
    return rows


def extend(model, sc, target=-4.):
    """Try the prior coordination slope, then adjacent SC increments."""
    nodes = [(0., np.zeros(3)), (-.5, sc)]
    slope_step = sc * .5
    increments = sorted((np.asarray(row) for row in product((-.25, 0., .25), repeat=3)),
                        key=lambda row: np.linalg.norm(row - slope_step))
    rejected = []
    for angle in np.arange(-.75, target - .01, -.25):
        start = nodes[-1]
        rejected = []
        for increment in increments:
            end = (float(angle), start[1] + increment)
            path = coordinated_segment(model, start, end)
            if path[-1]["hits"] or path[-1]["limit"]:
                rejected.append(dict(SC_rotvec_deg=end[1].tolist(), first_failure=path[-1]))
            else:
                nodes.append(end)
                break
        else:
            break
    ordered = list(zip(nodes, nodes[1:]))
    outbound = [coordinated_segment(model, a, b) for a, b in ordered]
    inbound = [coordinated_segment(model, b, a) for a, b in reversed(ordered)]
    return dict(last_angle_deg=nodes[-1][0], last_SC_rotvec_deg=nodes[-1][1].tolist(),
                outbound=outbound, inbound=inbound, rejected_next_step=rejected,
                outbound_and_inbound_clear=all(not segment[-1]["hits"] and not segment[-1]["limit"]
                                               for segment in (*outbound, *inbound)))


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1504_negative_rib_clearance.json"
    relief = read(source)
    coupled = load_reference(relief["source"])
    closure = load_reference(coupled["closure_audit"])
    design = load_reference(closure["source"])
    pose = load_reference(design["source"])
    access = load_reference(pose["source"])
    entry = load_reference(load_reference(access["source"])["source"])
    blind = load_reference(load_reference(entry["source"])["source"])
    split = load_reference(load_reference(load_reference(blind["parent"])["split_audit"])["split"])
    side = load_reference(load_reference(split["parent"])["parent"])
    removable = load_reference(side["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    sources = {key: meta for key, meta in assembly["sources"].items()
               if key not in entry["retired_keys"] and not (key.startswith("AC:") and key.endswith("_right"))}
    sources.update({key: meta for key, meta in entry["sources"].items() if key in sources})
    meshes = design["selected"]["meshes"]
    sources.update({"bone:scapula": meshes["bone:scapula"], "AC:AC_cap_right": meshes["cap_right"],
                    "AC:AC_screw_right": meshes["screw_right"],
                    "AC:AC_bottom_washer_left": access["left_washer"]["mesh"]})
    local = next(case["local_relief_meshes"] for case in relief["cases"] if case["target_deg"] == -.25)
    origins = load_reference(coupled["pose_reference"])["neutral_origins_m"]
    owners = {key: assembly["owners"].get(key, "scapula") for key in sources}
    cases = {}
    for name, parts in (("short_right_screw", sources), ("FJ1504_local_relief", {**sources, **local})):
        model = ShoulderAssembly(parts, owners, origins)
        zero = np.zeros(3)
        direct = trace(model, np.linspace(0., -1., 51), zero)
        safe = [row["scapular_coronal_deg"] for row in direct if not row["hits"] and not row["limit"]]
        reverse = trace(model, reversed(safe), zero)
        positive_path = coupled["cases"]["short_right_screw"]["directions"]["1"]["path"]
        positive = []
        for row in positive_path:
            hits, limit = model.evaluate(row["scapular_coronal_deg"], np.asarray(row["SC_rotvec_deg"]))
            if hits or limit:
                positive.append(dict(angle_deg=row["scapular_coronal_deg"], hits=hits, limit=limit))
        cases[name] = dict(parts=len(parts), negative_forward=direct, negative_reverse=reverse,
                           last_clear_deg=safe[-1] if safe else None,
                           previous_positive_path_failures=positive)
        print(name, cases[name]["last_clear_deg"], direct[-1]["scapular_coronal_deg"],
              direct[-1]["hits"], "positive path failures", len(positive), flush=True)
        if name == "FJ1504_local_relief":
            probes = []
            for vector in product((-.5, 0., .5), repeat=3):
                sc = np.array(vector)
                hits, limit = model.evaluate(-.5, sc)
                probe = dict(SC_rotvec_deg=sc.tolist(), hits=hits, limit=limit)
                if not hits and not limit:
                    for label, join in (("neutral", 0.), ("direct_limit", -.28)):
                        start = model.record(join, zero)
                        path = [dict(**start, hits=[], limit=None)]
                        for fraction in np.linspace(0., 1., 31)[1:]:
                            angle = join + fraction * (-.5 - join)
                            offset = sc * fraction
                            contact, bound = model.evaluate(angle, offset)
                            path.append(dict(**model.record(angle, offset), hits=contact, limit=bound))
                            if contact or bound:
                                break
                        probe[label + "_connection"] = path
                    if not path[-1]["hits"] and not path[-1]["limit"]:
                        initial = probe["neutral_connection"][-1]
                        if not initial["hits"] and not initial["limit"]:
                            probe["extension"] = extend(model, sc)
                            print("extended", probe["SC_rotvec_deg"],
                                  probe["extension"]["last_angle_deg"],
                                  probe["extension"]["outbound_and_inbound_clear"], flush=True)
                probes.append(probe)
            cases[name]["negative_half_degree_endpoint_probes"] = probes
            print("clear -0.5 endpoint probes", sum(not p["hits"] and not p["limit"] for p in probes),
                  flush=True)
    save(folder / "FJ1504_negative_relief_integration.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_negative_relief_integration.py"), source=provenance(source),
        revision_meshes=local, cases=cases, scope="Rigid assembly; SC, GH and elbow fixed at neutral "
        "on the negative trace. Cross-body overlap limit 0.001 mm3; 0.02 degree steps. "
        "Positive replay uses the previous reference samples. Flexible contact and loading are excluded.",
    ))


if __name__ == "__main__":
    main()
