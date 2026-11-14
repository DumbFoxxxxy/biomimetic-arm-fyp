"""Screen a separate proximal pad and belly on the revised FJ1500 route."""

import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import A, export, load_component, rod, save
from study_elbow_series_package import bounded, overlaps
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "s1_local_muscle_load/FJ1500_cap_relief"
    source = folder / "FJ1500_side_entry_relief_screen.json"
    relief = read(source)
    side = load_reference(relief["parent"])
    sleeve = load_reference(side["source"])
    feed = load_reference(sleeve["integrated_sleeve"])
    assembly = load_reference(feed["mount"])
    path = np.array(load_reference(feed["route"])["exported"]["trial"]["new_path_source_mm"])
    tpu = load_component(feed["selected"]["recessed_TPU"])
    cap = load_component(sleeve["selected"]["cap"])
    origin = path[2]
    axis = A.D.unit(path[-1]-origin)
    distal_edge = max(row["contact_station_range_mm"][1] for row in relief["checks"])
    rigid = {k:load_component(v) for k,v in assembly["sources"].items()
             if k != "cap:FJ1500/origin" and not k.startswith(("origin:FJ1504/","origin:FJ1500/"))}
    directions = [np.array(q["direction_source"]) for q in side["first_step_directions"][::4]]
    directions.extend((axis,-axis,np.array(feed["bore_axis_source"]),-np.array(feed["bore_axis_source"])))
    wire = A.D.outer_envelope(path,.30)
    stations = (path-origin)@axis

    def components(solid):
        return sorted((float(part.volume()) for part in solid.decompose() if part.volume()>1e-5),reverse=True)

    def approaches(moving, fixed):
        tested = []
        for direction in directions:
            worst = 0.
            for d in (0.25, .5, 1., 2., 4., 8.):
                hit = overlaps(moving.translate(d*direction),fixed)
                worst = max(worst,sum(hit.values()))
            tested.append((worst,direction))
        tested.sort(key=lambda row:row[0])
        dense = []
        for _, direction in tested[:3]:
            positions = []
            for distance in np.linspace(0.,8.,81):
                hits = overlaps(moving.translate(distance*direction),fixed)
                if hits:
                    positions.append(dict(distance_mm=float(distance),hits_mm3=hits))
            dense.append(dict(direction_source=direction.tolist(),hits=positions))
        return dict(directions=len(tested),best_coarse_overlap_mm3=tested[0][0],fine=dense)

    rows = []
    for split in (6.,10.,14.,18.,22.):
        level = float(origin@axis+split)
        pad = tpu.trim_by_plane(-axis,-level)
        belly = tpu.trim_by_plane(axis,level)
        reconstructed = pad+belly
        error = abs(float(reconstructed.volume()-tpu.volume()))
        if error>.001:
            raise ValueError("Split does not reconstruct source volume")
        indices = np.flatnonzero((stations[:-1]-split)*(stations[1:]-split)<0.)
        crossings = []
        for i in indices:
            fraction = (split-stations[i])/(stations[i+1]-stations[i])
            point = path[i]+fraction*(path[i+1]-path[i])
            crossings.append(point.tolist())
        row = dict(split_station_mm=split,channel_crossings_source_mm=crossings,
                   pad_component_volumes_mm3=components(pad),belly_component_volumes_mm3=components(belly),
                   source_volume_reconstruction_error_mm3=error,wire_overlap_mm3=float((wire^reconstructed).volume()))
        if len(row["pad_component_volumes_mm3"])==len(crossings)==1:
            row["pad_first"] = approaches(pad,bounded(rigid))
        if len(row["belly_component_volumes_mm3"])==len(crossings)==1:
            row["belly_without_pad"] = approaches(belly,bounded({**rigid,"cap":cap}))
            if any(not q["hits"] for q in row["belly_without_pad"]["fine"]):
                centre = np.array(crossings[0])
                reread, meta = export(belly,folder/f"FJ1500_split_belly_{split:g}mm.obj",origin)
                direction = np.array(next(q["direction_source"] for q in row["belly_without_pad"]["fine"] if not q["hits"]))
                targets = bounded({**rigid,"cap":cap})
                if any(overlaps(reread.translate(d*direction),targets) for d in np.linspace(0.,8.,81)):
                    raise ValueError("The belly approach did not reproduce after export")
                row["belly_mesh"] = meta
                # Probe the retained first millimetre for material to carry a
                # connector around the unchanged wire station.
                ring = rod(centre,centre+axis,3.)-rod(centre-.1*axis,centre+1.1*axis,.86)
                row["connector_stock_in_first_mm_mm3"] = float((ring^reread).volume())
        if "pad_first" in row and "belly_without_pad" in row:
            row["belly_second"] = approaches(belly,bounded({**rigid,"cap":cap,"proximal_pad":pad}))
        rows.append(row)
        print(split,len(row["pad_component_volumes_mm3"]),len(row["belly_component_volumes_mm3"]),len(crossings),
              [(key,row[key]["best_coarse_overlap_mm3"],min(len(q["hits"]) for q in row[key]["fine"]))
               for key in ("pad_first","belly_without_pad","belly_second") if key in row],flush=True)
    save(folder/"FJ1500_split_origin_screen.json",dict(
        generator=provenance(ROOT/"tools/study_s1_infraspinatus_split_origin.py"), parent=provenance(source),
        observed_contact_distal_station_mm=distal_edge, cases=rows,
        scope="Unfastened planar split at the existing lumen. Nominal rigid insertion geometry."
              " The proximal pad is screened without the cap. Exported belly approaches omit the pad and sleeve."
              " Cap installation, fasteners, a load-bearing pad-to-belly connection and TPU flexibility remain to be designed."))


if __name__ == "__main__":
    main()
