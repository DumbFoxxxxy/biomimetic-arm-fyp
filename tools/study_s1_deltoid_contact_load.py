"""Assign FJ1468 guide reactions to actual surfaces and screen wall compliance."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, load_component, save
from channel_first_bone_contact import BoneContactSurface
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1468_taut_wire_neutral.json"
    taut = read(source)
    motion = load_reference(taut["source"])
    audit = load_reference(motion["wall_audit"])
    transition = load_reference(audit["source"])
    proximal = load_reference(load_reference(transition["source"])["proximal_mount"])["sources"]
    solids = {
        "TPU": load_component(motion["merged_TPU"]),
        "cap": load_component(proximal["cap:FJ1468/origin"]),
    }
    best = max(taut["cases"], key=lambda row: row["node_count"])
    if not best["solver_success"] or best["full_wire_hits_mm3"] or not best["reaction_samples"]:
        raise ValueError("The fine taut-wire case has no valid reaction samples")
    centres = np.asarray([row["point_source_mm"] for row in best["reaction_samples"]])
    lo, hi = centres.min(axis=0) - 3., centres.max(axis=0) + 3.
    region = A.I._manifold().Manifold.cube(hi - lo).translate(lo)
    surfaces = {name: BoneContactSurface(part ^ region, band_mm=1., cell_mm=1.)
                for name, part in solids.items()}
    distances = {name: surface.query(centres) for name, surface in surfaces.items()}
    rows = []
    for i, sample in enumerate(best["reaction_samples"]):
        gap = {name: float(values[0][i]) for name, values in distances.items()}
        contact_owners = [name for name, distance in gap.items() if abs(distance - .3) <= .01]
        owner = contact_owners[0] if len(contact_owners) == 1 else "unresolved"
        normal = distances[owner][1][i] if owner != "unresolved" else np.zeros(3)
        reaction = sample["normal_reaction_N_at_10N"]
        rows.append(dict(**sample, surface_distance_mm=gap, owner=owner,
                         wall_point_source_mm=(centres[i] - gap[owner] * normal).tolist()
                         if owner != "unresolved" else None,
                         force_on_wall_N_at_10N=(-reaction * normal).tolist()
                         if owner != "unresolved" else None))
    summaries = {}
    for name in (*solids, "unresolved"):
        matching = [row for row in rows if row["owner"] == name]
        summaries[name] = dict(samples=len(matching),
                               scalar_reaction_N_at_10N=sum(row["normal_reaction_N_at_10N"] for row in matching),
                               resultant_N_at_10N=np.sum([row["force_on_wall_N_at_10N"] for row in matching
                                                          if row["force_on_wall_N_at_10N"] is not None],
                                                         axis=0).tolist() if matching and name != "unresolved"
                               else None)
    # Sum vectors first. This is an equivalent one-DOF stiffness of the loaded
    # TPU region as a whole, independent of contact sample count.
    force = float(np.linalg.norm(summaries["TPU"]["resultant_N_at_10N"] or [0., 0., 0.]))
    compliance = [dict(line_tension_N=tension, patch_stiffness_N_per_mm=k,
                       predicted_patch_displacement_mm=force * tension / (10. * k))
                  for tension in (5., 10., 20.) for k in (5., 10., 20., 40.)]
    save(folder / "FJ1468_contact_load_screen.json", dict(
        generator=provenance(ROOT / "tools/study_s1_deltoid_contact_load.py"),
        source=provenance(source), rows=rows, by_owner=summaries, compliance=compliance,
        scope="Neutral, fixed-wall, frictionless contact reactions from the fine discrete cable solve. "
              "Ownership is inferred by proximity to the cable radius. The compliance table uses an "
              "unmeasured equivalent patch stiffness with reactions held fixed; it is not a coupled "
              "solution, local indentation estimate, or a prediction of available wire stroke.",
    ))
    print(summaries, flush=True)
    print("10 N displacement cases", [row for row in compliance if row["line_tension_N"] == 10.])


if __name__ == "__main__":
    main()
