"""Fit a removable flanged guide in the FJ1500 scapular cap."""

from itertools import product

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_feed_sleeve_screen.json"
    feed = read(source)
    assembly = load_reference(feed["mount"])
    route = load_reference(feed["route"])["exported"]["trial"]
    path = np.array(route["new_path_source_mm"])
    entry = np.array(feed["start_source_mm"])
    axis = np.array(feed["bore_axis_source"])
    end = feed["span_mm"]-.1
    cap_key = "cap:FJ1500/origin"
    cap = load_component(feed["selected"]["mesh"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    fixed = bounded({key: load_component(src) for key, src in assembly["sources"].items() if key != cap_key})
    wire = A.D.outer_envelope(path,.30)

    def axial_rod(lo, hi, radius):
        return rod(entry+lo*axis,entry+hi*axis,radius,64)

    trials = []
    selected = None
    for radius, seat, flange_radius in product((1.35,1.41),(0.,.5,1.,2.),(1.85,2.15,2.45)):
        body = axial_rod(-.8,end,radius)
        flange = axial_rod(seat-.8,seat,flange_radius)
        solid = ((body+flange)-axial_rod(-1.,end+.1,.86)).as_original().simplify(1e-6)
        pocket = axial_rod(-9.,end+.1,radius+.05)+axial_rod(-9.,seat,flange_radius+.05)
        revised = (cap-pocket).as_original().simplify(1e-6)
        components = [sum(q.volume()>1e-5 for q in piece.decompose()) for piece in (revised,solid)]
        stationary = {**fixed,**bounded({"cap":revised,"TPU":tpu})}
        collisions = overlaps(solid,stationary)
        seat_witness = float((solid.translate(.05*axis)^revised).volume())
        row = dict(sleeve_outer_radius_mm=radius,flange_radius_mm=flange_radius,seat_station_mm=seat,
                   components=components,installed_hits_mm3=collisions, seat_witness_mm3=seat_witness,
                   nominal_wire_overlap_mm3=float((wire^(revised+solid)).volume()),
                   removed_cap_volume_mm3=float((cap-revised).volume()))
        trials.append(row)
        if components != [1,1] or collisions or seat_witness<.05 or row["nominal_wire_overlap_mm3"]>.001:
            continue
        row["approach_hits"] = []
        for distance in np.linspace(8.,0.,81):
            hits = overlaps(solid.translate(-distance*axis),stationary)
            if hits:
                row["approach_hits"].append(dict(offset_mm=float(distance),hits_mm3=hits))
        if row["approach_hits"] or selected is not None:
            continue
        reread, cap_source = export(revised,folder/"FJ1500_removable_feed_cap.obj",entry)
        sleeve, sleeve_source = export(solid,folder/"FJ1500_removable_feed_sleeve.obj",entry)
        targets = {**fixed,**bounded({"cap":reread,"TPU":tpu})}
        if (any(overlaps(sleeve.translate(-distance*axis),targets) for distance in np.linspace(8.,0.,81))
                or (sleeve.translate(.05*axis)^reread).volume()<.05 or (wire^(reread+sleeve)).volume()>.001):
            raise ValueError("Fit did not reproduce after mesh export")
        selected = dict(**row,cap=cap_source,sleeve=sleeve_source)
    tpu_installation = []
    if selected is not None:
        assembly_file = folder/"FJ1500_feed_installation_screen.json"
        assembly_screen = read(assembly_file)
        if assembly_screen["source"] != provenance(source):
            raise ValueError("Installation directions use a different feed model")
        targets = {**{key:part for key,part in fixed.items() if not key.startswith("origin:FJ1504/")},
                   **bounded({"cap":load_component(selected["cap"])})}
        seated = overlaps(tpu,targets)
        directions = {row["direction"]:row["vector_source"] for row in assembly_screen["results"]
                      if row["sequence"] == "mount_TPU_after_cap"}
        for label, vector in directions.items():
            hits = []
            for distance in np.linspace(8.,0.,41):
                current = overlaps(tpu.translate(distance*np.array(vector)),targets)
                if current:
                    hits.append(dict(offset_mm=float(distance),hits_mm3=current))
            tpu_installation.append(dict(direction=label,vector_source=vector,seated_hits_mm3=seated,hits=hits))
    save(folder/"FJ1500_removable_feed_screen.json",dict(
        generator=provenance(ROOT/"tools/build_s1_infraspinatus_removable_feed.py"), integrated_sleeve=provenance(source),
        bore_axis_source=axis.tolist(), insertion_direction_source=(-axis).tolist(), trials=trials,selected=selected,
        tpu_approaches_with_sleeve_removed=tpu_installation,
        scope="Flange transfers axial compression to the counterbore. Axial stop witness is a solid translation,"
              " not a load test. Sleeve insertion assumes TPU already seated. Reverse retention, complete entry"
              " beyond the sampled 8 mm, TPU insertion and capacity remain unverified."))
    print('Trials',len(trials),'selected',bool(selected),flush=True)


if __name__ == "__main__":
    main()
