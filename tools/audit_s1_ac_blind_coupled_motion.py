"""Compare SC/AC coordination before and after the short AC screw revision."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from study_s1_coupled_scapula import ShoulderAssembly, follow
from study_s1_local_force_geometry import load_reference


def load_assembly(folder):
    audit_file = folder / "FJ1504_AC_blind_closure_audit.json"
    audit = read(audit_file)
    design = load_reference(audit["source"])
    pose = load_reference(design["source"])
    access = load_reference(pose["source"])
    entry = load_reference(load_reference(access["source"])["source"])
    cap = load_reference(entry["source"])
    blind = load_reference(cap["source"])
    split = load_reference(load_reference(load_reference(blind["parent"])["split_audit"])["split"])
    removable = load_reference(load_reference(load_reference(split["parent"])["parent"])["source"])
    assembly = load_reference(load_reference(removable["integrated_sleeve"])["mount"])
    sources = {key: meta for key, meta in assembly["sources"].items()
               if key not in entry["retired_keys"]}
    sources.update({key: meta for key, meta in entry["sources"].items() if key in sources})
    revised = {key: meta for key, meta in sources.items()
               if not (key.startswith("AC:") and key.endswith("_right"))}
    meshes = design["selected"]["meshes"]
    revised.update({"bone:scapula": meshes["bone:scapula"], "AC:AC_cap_right": meshes["cap_right"],
                    "AC:AC_screw_right": meshes["screw_right"],
                    "AC:AC_bottom_washer_left": access["left_washer"]["mesh"]})
    return assembly, sources, revised


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    audit_file = folder / "FJ1504_AC_blind_closure_audit.json"
    assembly, sources, revised = load_assembly(folder)
    passive = DEST / "elbow_passive_model/model_check.json"
    origins = read(passive)["neutral_origins_m"]
    outcomes = {}
    for name, parts in (("original_AC", sources), ("short_right_screw", revised)):
        owners = {key: assembly["owners"].get(key, "scapula") for key in parts}
        model = ShoulderAssembly(parts, owners, origins)
        neutral, _ = model.evaluate(0., np.zeros(3))
        result = dict(rigid_part_count=len(parts), neutral_hits=neutral, directions={})
        if not neutral:
            for sign in (-1, 1):
                path = follow(model, 20, sign, True)
                reverse = []
                for row in reversed(path["path"]):
                    angle = row["scapular_coronal_deg"]
                    sc = np.asarray(row["SC_rotvec_deg"])
                    hits, limit = model.evaluate(angle, sc)
                    if hits or limit:
                        reverse.append(dict(angle_deg=angle, hits=hits, limit=limit))
                path["reverse_replay_hits"] = reverse
                result["directions"][str(sign)] = path
                print(name, sign, path["path"][-1]["scapular_coronal_deg"], len(reverse), flush=True)
            if name == "short_right_screw" and len(result["directions"]["-1"]["path"]) == 1:
                service = []
                zero = np.zeros(3)
                for vector in sorted(product((-.5, 0., .5), repeat=3), key=np.linalg.norm):
                    sc = np.asarray(vector)
                    preposition, blocked = model.transition(0., 0., zero, sc)
                    row = dict(preposition_deg=sc.tolist(), preposition_failure=blocked)
                    if blocked is None:
                        negative, failure = model.transition(0., -.5, sc, sc)
                        row.update(preposition_path=preposition, negative_path=negative, negative_failure=failure)
                    service.append(row)
                result["SC_preposition_probes"] = service
                print("SC preposition", sum(row["preposition_failure"] is None for row in service),
                      "negative half-degree paths",
                      sum(row.get("negative_failure", True) is None for row in service), flush=True)
        else:
            print(name, "neutral blocked", neutral[:6], flush=True)
        outcomes[name] = result
    save(folder / "FJ1504_AC_blind_coupled_motion.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_ac_blind_coupled_motion.py"),
        closure_audit=provenance(audit_file), pose_reference=provenance(passive),
        cases=outcomes, scope="Ideal SC/AC rigid geometry with GH and elbow held at neutral, "
        "and prescribed coronal scapular rotation. No flexible parts or driving forces. "
        "Discrete greedy path search and reverse replay; not a global feasibility proof.",
    ))


if __name__ == "__main__":
    main()
