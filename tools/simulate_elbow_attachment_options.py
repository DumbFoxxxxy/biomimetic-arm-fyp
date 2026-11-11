"""Compare shortlisted force application points with the existing stop model."""

from pathlib import Path

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import HEADS, JOINTS, measure, simulate, write_csv


def main():
    output = DEST.parent / "elbow_local_muscle_load"
    options_file = output / "elbow_local_attachment_options.json"
    stop_file = output / "elbow_local_stop_response.json"
    options, stop = read(options_file), read(stop_file)
    model_path = ROOT / stop["local_model"]["file"]
    if provenance(model_path) != stop["local_model"]:
        raise RuntimeError("Local stop model fingerprint changed")
    xml = model_path.read_text(encoding="utf-8")
    cases = []
    for arrangement in ("original", "shortlist"):
        model = mujoco.MjModel.from_xml_string(xml)
        data = mujoco.MjData(model)
        mujoco.mj_forward(model, data)
        addresses = [model.joint(name).qposadr.item() for name in JOINTS]
        sites = [[model.site(f"{head}_{end}").id for end in ("origin", "distal")] for head in HEADS]
        positions = []
        for head, (_, distal) in zip(HEADS, sites, strict=True):
            if arrangement == "shortlist":
                key = ",".join(map(str, options["local_probe_shortlist"][head]))
                candidate = options["tested_options"][head][key]
                position = np.asarray(candidate["point_source_mm"]) * .001
                body = data.body(int(model.site_bodyid[distal]))
                model.site_pos[distal] = body.xmat.reshape(3, 3).T @ (position - body.xpos)
                positions.append(position)
            else:
                positions.append(data.site_xpos[distal].copy())
        mujoco.mj_forward(model, data)
        np.testing.assert_allclose(data.site_xpos[[row[1] for row in sites]], positions, atol=1e-9)
        neutral, _ = measure(model, data, sites)
        for ext_stiffness in (200., 50.):
            trace, metrics = simulate(model, sites, neutral, addresses, (200., ext_stiffness), 25., 0, .0025)
            path = output / f"elbow_attachment_{arrangement}_ext{ext_stiffness:g}.csv"
            write_csv(path, ["time_s", *(f"{joint}_deg" for joint in JOINTS),
                             *(f"{head}_length_change_mm" for head in HEADS),
                             *(f"{head}_drive_N" for head in HEADS),
                             *(f"{head}_TPU_restoring_N" for head in HEADS)], trace)
            cases.append(dict(**metrics, arrangement=arrangement, trace=provenance(path)))
            print(arrangement, ext_stiffness, "peak", round(metrics["elbow_max_deg"], 4),
                  "peak cap rotations", metrics["cap_peak_abs_deg"],
                  "extensor length gain", metrics["muscle_max_length_change_mm"][1], flush=True)
    save(output / "elbow_attachment_response.json", dict(
        generator=provenance(Path(__file__).resolve()), stop=provenance(stop_file),
        options=provenance(options_file), local_model=provenance(model_path),
        stop_reference=stop["stop_reference"], cases=cases,
        scope="Force-only comparison. Candidate bore coordinates are hypothetical, attached to the "
        "existing moving cap frames. Model mass and pivots remain original; no candidate connector or "
        "rigid/TPU collision geometry is included. Springs start at each neutral span, with linear "
        "tension and compression stiffness equal within each muscle. Caps and elbow move freely."))


if __name__ == "__main__":
    main()
