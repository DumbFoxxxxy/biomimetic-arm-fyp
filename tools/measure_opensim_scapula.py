"""Measure constrained OpenSim shoulder path geometry and moment arms."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import opensim as osim


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "02_人体骨骼晶格肌肉方案/输出/whole_arm_first_pass/scapular_control"
GROUPS = ("SerratusAnterior_", "Trapezius", "Rhomboideus_", "LevatorScapulae", "PectoralisMinor")
COORDINATES = ("scapula_upward_rot", "scapula_abduction", "scapula_elevation", "scapula_winging")


def xyz(value):
    return [1000. * value.get(index) for index in range(3)]


def sample(model, coordinate, baseline, offset, muscles):
    # Make the key coordinate exact while assembling AC closure. Unlock it
    # again before requesting the native muscle moment arm.
    request = baseline + math.radians(offset)
    coordinate.setDefaultValue(request)
    coordinate.setDefaultLocked(True)
    state = model.initSystem()
    model.assemble(state)
    coordinate.setLocked(state, False)
    model.realizePosition(state)
    achieved = coordinate.getValue(state)
    if abs(request - achieved) > math.radians(.1):
        raise ValueError(f"Assembled angle differs from requested angle {offset}")
    records = {}
    for muscle in muscles:
        path = muscle.getGeometryPath()
        current = path.getCurrentPath(state)
        points = []
        for index in range(current.getSize()):
            point = current.get(index)
            points.append(dict(name=point.getName(), frame=point.getParentFrame().getName(),
                               type=point.getConcreteClassName(), body_local_mm=xyz(point.getLocation(state)),
                               ground_mm=xyz(point.getLocationInGround(state))))
        records[muscle.getName()] = dict(length_mm=1000. * path.getLength(state),
                                         native_moment_arm_mm=1000. * muscle.computeMomentArm(state, coordinate),
                                         points=points)
    values = model.getCoordinateSet()
    other = {name: math.degrees(values.get(name).getValue(state))
             for name in ("clav_prot", "clav_elev", *COORDINATES)}
    return dict(requested_offset_deg=offset, offset_deg=math.degrees(achieved - baseline),
                coordinate_deg=math.degrees(achieved), angles_deg=other, muscles=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "OpenSimDownloads/TSM_subject_noWeight.osim")
    args = parser.parse_args()
    path = args.model.resolve()
    geometry = path.parent / "Geometry"
    if geometry.is_dir():
        osim.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry))
    model = osim.Model(str(path))
    state = model.initSystem()
    model.assemble(state)
    coordinate = model.getCoordinateSet().get("scapula_upward_rot")
    baseline = coordinate.getValue(state)
    muscles = model.getMuscles()
    selected = [muscles.get(i) for i in range(muscles.getSize()) if muscles.get(i).getName().startswith(GROUPS)]
    outbound = [sample(model, coordinate, baseline, i / 2, selected) for i in range(21)]
    returning = [sample(model, coordinate, baseline, i / 2, selected) for i in range(20, -1, -1)]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    result = dict(source=str(path.relative_to(ROOT)), source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  version=osim.GetVersion(), coordinate=coordinate.getName(), model_neutral_deg=math.degrees(baseline),
                  step_deg=.5, outbound=outbound, returning=returning,
                  scope="Original joint and AC constraint; the selected coordinate is prescribed. "
                        "Moment arms are native OpenSim; slopes of the assembled sweep are reported separately.")
    (OUTPUT / "opensim_scapula.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                                                encoding="utf-8", newline="\n")
    with (OUTPUT / "opensim_scapula.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("direction", "offset_deg", "muscle", "length_mm", "native_moment_arm_mm",
                         "negative_length_gradient_mm"))
        for direction, rows in (("outbound", outbound), ("return", returning)):
            for index, row in enumerate(rows):
                left, right = rows[max(0, index - 1)], rows[min(len(rows) - 1, index + 1)]
                delta = math.radians(right["offset_deg"] - left["offset_deg"])
                for name, item in row["muscles"].items():
                    gradient = (right["muscles"][name]["length_mm"] - left["muscles"][name]["length_mm"]) / delta
                    writer.writerow((direction, row["offset_deg"], name, item["length_mm"],
                                     item["native_moment_arm_mm"], -gradient))
    print(f"OpenSim: {len(selected)} bundles, {len(outbound)} outward and {len(returning)} return samples.")


if __name__ == "__main__":
    main()
