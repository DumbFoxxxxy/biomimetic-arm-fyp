"""Export the neutral OpenSim comparison in the first-pass web mesh frame."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "02_人体骨骼晶格肌肉方案/输出/whole_arm_first_pass"
VIEWER = ROOT / "docs/viewer"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = FOLDER / "opensim_comparison.json"
    report = read(source)
    manifest = read(VIEWER / "first-pass.json")
    assembly = report["engineering"]
    if (assembly["sha256"] != sha256(ROOT / assembly["file"])
            or assembly["sha256"] != manifest["assembly_sha256"]
            or assembly["file"] != manifest["assembly_source"]):
        raise ValueError("OpenSim comparison and web geometry use different assemblies")
    guide_ref = report["guide_points"]
    if guide_ref["sha256"] != sha256(ROOT / guide_ref["file"]):
        raise ValueError("OpenSim comparison uses an older guide-point record")
    guides = read(ROOT / guide_ref["file"])["data"]
    origins = read(FOLDER / "whole_arm_shared_scene_check.json")["neutral_body_origins_mm"]
    pivot = manifest["source_origin_mm"]

    def to_web(point, owner=None):
        origin = origins[owner] if owner else [0., 0., 0.]
        return [round((x + offset - center) * .001, 8) for x, offset, center in zip(point, origin, pivot, strict=True)]

    routes = {}
    for key, row in report["routes"].items():
        guide = guides[key]
        loaded = "centreline_source_mm" not in guide
        line = guide["loaded_neutral_tendon_tips_world_mm"] if loaded else guide["centreline_source_mm"]
        paths, candidates = [], []
        for bundle in row["anatomical_paths"]:
            segment = []
            for point in bundle["points"]:
                if "surface_local_mm" in point:
                    candidates.append(dict(role=point["role"], xyz=to_web(point["surface_local_mm"], point["owner"])))
                if point["type"] == "ConditionalPathPoint":
                    low, high = map(float, point["active_range"].split())
                    if not low <= 0. <= high:
                        continue
                if "mapped_local_mm" in point:
                    segment.append(to_web(point["mapped_local_mm"], point["owner"]))
                else:
                    if len(segment) > 1:
                        paths.append(segment)
                    segment = []
            if len(segment) > 1:
                paths.append(segment)
        routes[key] = dict(label=row["label"], zone=row["anatomical_zone"],
                           source_model="Arm26" if key in ("FJ1486", "FJ1477") else "Seth / Belli",
                           bundles=[bundle["muscle"] for bundle in row["anatomical_paths"]],
                           reference=paths, candidates=candidates, engineering=[to_web(point) for point in line],
                           engineering_ports=[to_web(point, owner) for point, owner in
                                              zip(row["engineering_ports_local_mm"], row["engineering_owners"], strict=True)],
                           support=row["engineering_support"], loaded_chord_10n=loaded)
    summary = dict(route_count=len(routes), bundle_count=sum(len(route["bundles"]) for route in routes.values()))
    data = dict(schema_version=1, units="m", assembly_sha256=assembly["sha256"], mesh_sha256=manifest["mesh_sha256"],
                report_sha256=sha256(source), neutral=True, summary=summary, routes=routes)
    (VIEWER / "opensim-overlay.json").write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n",
                                                encoding="utf-8", newline="\n")
    print(f'Exported {summary["route_count"]} routes and {summary["bundle_count"]} anatomical bundles')


if __name__ == "__main__":
    main()
