"""Check three continuous SC/AC and GH flexion paths on the revised rigid assembly."""

from itertools import combinations

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_series_elbow_integration import MOBILE, compose, find_hits
from audit_s1_girdle_tendon_travel import place
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_local_force_geometry import load_reference


def audit_path(path, evaluate):
    cleared = []
    first = None
    for pose in path:
        hits = evaluate(pose)
        if hits:
            low = np.array(cleared[-1]) if cleared else None
            high = np.array(pose)
            if low is not None:
                while np.max(np.abs(high - low)) > .1:
                    middle = (low + high) / 2.
                    middle_hits = evaluate(middle)
                    if middle_hits:
                        high, hits = middle, middle_hits
                    else:
                        low = middle
                cleared.append(low.tolist())
            first = dict(pose_deg=high.tolist(), hits=hits)
            break
        cleared.append(list(pose))
    return dict(clear_path_deg=cleared, first_hit=first,
                reverse_replay_hits=[dict(pose_deg=pose, hits=hits) for pose in reversed(cleared)
                                     if (hits := evaluate(pose))])


def main():
    local = DEST.parent / "s1_local_muscle_load"
    folder = local / "FJ1500_cap_relief"
    path_ref = provenance(folder / "FJ1504_AC_blind_coupled_motion.json")
    span_ref = provenance(local / "FJ1468_girdle_flexion_span.json")
    rigid = load_reference(path_ref)
    spans = load_reference(span_ref)
    if spans["rigid_path"] != path_ref or rigid["closure_audit"] != provenance(folder / "FJ1504_AC_blind_closure_audit.json"):
        raise ValueError("Shoulder geometry inputs changed")
    endpoints = load_reference(spans["local_endpoints"])
    origins = load_reference(rigid["pose_reference"])["neutral_origins_m"]
    assembly, _, sources = load_assembly(folder)
    owners = {key: assembly["owners"].get(key, "scapula") for key in sources}
    model = ShoulderAssembly(sources, owners, origins)
    candidate = rigid["cases"]["short_right_screw"]
    outward = candidate["directions"]["1"]
    if len(sources) != candidate["rigid_part_count"] or outward["reverse_replay_hits"] or candidate["neutral_hits"]:
        raise ValueError("Expected checked short-screw assembly")
    if not np.allclose(model.gh, endpoints["gh_centre_source_mm"], atol=.001, rtol=0.):
        raise ValueError("Mismatched GH reference")
    groups = [3 if owners[key] in MOBILE else group for key, group in zip(model.keys, model.groups, strict=True)]
    pairs = np.array([(i, j) for i, j in combinations(range(len(groups)), 2) if groups[i] != groups[j]])
    axis = np.asarray(endpoints["axes_source"][0])
    path = np.array([[row["scapular_coronal_deg"], *row["SC_rotvec_deg"]] for row in outward["path"]])
    stop = path[-1, 0]

    def evaluate(pose):
        scapula_deg, *sc, flexion_deg = pose
        clavicle, scapula, _ = model.kinematics(scapula_deg, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flexion_deg) * axis).as_matrix(), model.gh)
        frames = (None, clavicle, scapula, compose(scapula, gh))
        placed = [part.transform(frames[group]) if group else part
                  for part, group in zip(model.solids, groups, strict=True)]
        return find_hits(placed, model.keys, pairs)

    results = {}
    gh_samples = spans["gh_flexion_grid_deg"]
    if np.max(np.abs(np.diff(path, axis=0))) > 1. or np.max(np.diff(gh_samples)) > 1. or stop <= 0.:
        raise ValueError("Unexpected sampling grid")
    simultaneous = []
    initial = np.r_[path[0], 0.]
    simultaneous.append(initial)
    for row in path[1:]:
        target = np.r_[row, gh_samples[-1] * row[0] / stop]
        count = max(1, int(np.ceil(np.max(np.abs(target - initial)))))
        simultaneous.extend(initial + fraction * (target - initial) for fraction in np.linspace(0., 1., count + 1)[1:])
        initial = target
    for name, sequence in (
        ("gh_only", [np.r_[path[0], angle] for angle in gh_samples]),
        ("girdle_then_gh", [*[np.r_[row, 0.] for row in path], *[np.r_[path[-1], angle] for angle in gh_samples[1:]]]),
        ("simultaneous", simultaneous),
    ):
        results[name] = audit_path(sequence, evaluate)
        print(name, results[name]["clear_path_deg"][-1:], results[name]["first_hit"], flush=True)
    save(local / "FJ1468_girdle_flexion_clearance.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_girdle_flexion_clearance.py"),
        girdle_path=path_ref, span_study=span_ref, rigid_part_count=len(groups), cross_group_pairs=len(pairs),
        pose_columns=["scapular_coronal", "SC_x", "SC_y", "SC_z", "GH_flexion"],
        cases=results, scope="Prescribed sampled paths, rigid geometry only, elbow held neutral. "
                             "Regular parameter increments at most 1 degree, first hit bracket at most 0.1 degree. "
                             "Return retraces the clear samples. Flexible hardware and force response are excluded.",
    ))


if __name__ == "__main__":
    main()
