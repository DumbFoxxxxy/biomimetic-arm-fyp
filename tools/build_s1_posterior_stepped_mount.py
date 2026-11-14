"""Build a conforming FJ1467 blind mount on the current S1 scapula."""

from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from audit_s1_ac_blind_coupled_motion import load_assembly
from audit_s1_series_elbow_integration import MOBILE, compose
from build_channel_first_inventory import DEST, provenance
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps, union
from study_s1_coupled_scapula import ShoulderAssembly, affine_at
from study_s1_local_force_geometry import load_reference


def main():
    local = DEST.parent / "s1_local_muscle_load"
    stock_ref = provenance(local / "FJ1467_current_scapula_blind_axes.json")
    stock = load_reference(stock_ref)
    locations = load_reference(stock["route"])
    neighbor = load_reference(locations["neighbor_replay"])
    rigid = load_reference(neighbor["rigid_audit"])
    girdle = load_reference(rigid["girdle_path"])
    endpoints = load_reference(load_reference(rigid["span_study"])["local_endpoints"])
    assembly, _, sources = load_assembly(local / "FJ1500_cap_relief")
    bone_key = "bone:scapula"
    if sources[bone_key] != stock["scapula"]:
        raise ValueError("Scapula changed after axis screening")
    sources["elbow:FJ1486_bone_adapter"] = neighbor["adapter_mesh"]
    original = {key: load_component(ref) for key, ref in sources.items()}
    bone = original[bone_key]
    normal = np.asarray(stock["normal_source"])
    port = np.asarray(locations["source_ports"]["origin"]["inner_face_source_mm"])
    cable = np.asarray(locations["neutral_source_route_mm"])
    selected = stock["compact_pair_nearest_centre"]
    mounts = [selected[side] for side in ("first", "second")]
    mf = A.I._manifold()
    bridge = mf.Manifold.batch_hull([rod(point, point + 2. * normal, 3.2, 48)
                                   for point in (np.asarray(row["axis_source_mm"]) for row in mounts)])
    feet, bores, fasteners = [], [], {}
    for index, mount in enumerate(mounts):
        axis = np.asarray(mount["axis_source_mm"])
        surface = np.asarray(mount["bone_entry_mm"])
        feet.append(rod(surface - 3. * normal, axis + 1. * normal, 2.6, 48))
        bores.append(rod(surface - 3. * normal, axis + 2.1 * normal, 1.15, 48))
        shaft = rod(surface - 2.8 * normal, axis + 2.05 * normal, .95, 48)
        head = rod(axis + 2. * normal, axis + 3. * normal, 1.65, 48)
        fasteners[f"fastener_{index}"] = shaft + head
    added = union([bridge, *feet])
    relieved_bone = bone - union(bores)
    lumen = A.D.outer_envelope(cable[:12], .86)
    # Trim the underside directly against the recipient bone. Drill only
    # the fastener bores into bone; the wire lumen is formed in the support.
    support = ((added - bone) - union([lumen, *bores])).as_original().simplify(1e-7)
    parts = dict(support=support, **fasteners)
    candidates = {bone_key: relieved_bone, **parts}
    model = ShoulderAssembly({}, {}, load_reference(girdle["pose_reference"])["neutral_origins_m"])
    axis = np.asarray(endpoints["axes_source"][0])
    owners = {key: assembly["owners"].get(key, "scapula") for key in original}
    components = {key: [float(piece.volume()) for piece in solid.decompose()] for key, solid in candidates.items()}
    meshes = {}
    if all(sum(v > 1e-8 for v in vv) == 1 for vv in components.values()):
        for key, solid in candidates.items():
            name = key.replace(":", "_")
            _, meshes[key] = export(solid, local / f"FJ1467_stepped_{name}.obj", port)
        candidates = {key: load_component(meta) for key, meta in meshes.items()}
        relieved_bone = candidates[bone_key]
        support = candidates["support"]
        fasteners = {key: candidates[key] for key in fasteners}
        parts = dict(support=support, **fasteners)
    rows = []
    for sample in locations["cases"]:
        angle, *sc, flex = sample["pose_deg"]
        clavicle, scapula, _ = model.kinematics(angle, np.asarray(sc))
        gh = affine_at(Rotation.from_rotvec(np.deg2rad(flex) * axis).as_matrix(), model.gh)
        humerus = compose(scapula, gh)
        frames = {"clavicle": clavicle, "scapula": scapula, **dict.fromkeys(MOBILE, humerus)}
        targets = bounded({key: solid.transform(frames[owners[key]]) if owners[key] in frames else solid
                           for key, solid in {**original, bone_key: relieved_bone}.items()})
        hits = {key: hit for key, solid in parts.items()
                if (hit := overlaps(solid.transform(scapula), targets))}
        rows.append(dict(pose_deg=sample["pose_deg"], rigid_hits_mm3=hits))
    own_hits = {key: float((part ^ support).volume()) for key, part in fasteners.items()}
    local_wire = A.D.outer_envelope(cable[:12], .30)
    passage = overlaps(local_wire, bounded(candidates))
    save(local / "FJ1467_stepped_mount.json", dict(
        generator=provenance(Path(__file__)), stock_study=stock_ref, components_mm3=components,
        bone_removed_mm3=float((bone - relieved_bone).volume()), meshes=meshes,
        fastener_support_overlap_mm3=own_hits, neutral_proximal_wire_hits_mm3=passage, cases=rows,
        scope="Nominal rigid envelopes. Fasteners are shaft and head solids, without thread retention "
              "or preload. Support is conformed to the current scapula with two blind bores. "
              "Cable bone passage, TPU attachment, and installation tool access require separate checks.",
    ))
    print("Components:", {key: len(value) for key, value in components.items()},
          "blocked poses:", sum(bool(row["rigid_hits_mm3"]) for row in rows),
          "neutral wire:", passage, flush=True)


if __name__ == "__main__":
    main()
