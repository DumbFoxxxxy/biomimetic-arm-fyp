"""Repeat local equivalent muscle loading against the existing stop-face pair."""

import copy
from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import HEADS, JOINTS, build_local_model, measure, simulate, write_csv


def main():
    folder = DEST / "elbow_extension_stop"
    stop_check_path = folder / "model_check.json"
    stop_check = read(stop_check_path)
    passive_path = ROOT / stop_check["parent"]["file"]
    passive = read(passive_path)
    stop_source = ROOT / stop_check["model"]["file"]
    passive_source = ROOT / passive["model"]["file"]
    for path, expected in ((stop_source, stop_check["model"]), (passive_source, passive["model"]),
                           (passive_path, stop_check["parent"])):
        if provenance(path) != expected:
            raise RuntimeError(f"Source fingerprint changed: {path.name}")

    stop = ET.parse(stop_source)
    local = ET.fromstring(build_local_model(passive_source, passive))
    for name in ("forearm", "humerus"):
        body = local.find(f".//body[@name='{name}']")
        body.remove(body.find("inertial"))
        body.insert(0, copy.deepcopy(stop.find(f".//body[@name='{name}']/inertial")))
    for name, owner in (("fixed_extension_face", "humerus"), ("moving_extension_face", "forearm")):
        local.find(f".//body[@name='{owner}']").append(copy.deepcopy(stop.find(f".//geom[@name='{name}']")))
    local.find(".//joint[@name='elbow']").set("limited", "false")
    ET.SubElement(local, "contact").append(copy.deepcopy(stop.find("./contact/pair[@name='actual_extension_stop']")))
    option = local.find("option")
    option.attrib.update(stop.find("option").attrib)
    output = DEST.parent / "elbow_local_muscle_load"
    xml_path = output / "elbow_local_force_with_stop.xml"
    xml = ET.tostring(local, encoding="unicode")
    xml_path.write_text(xml, encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    if model.nv != 3 or model.npair != 1 or model.joint("elbow").limited.item():
        raise RuntimeError("Expected the three-joint local model with just the extension stop face pair")
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    site_ids = [[model.site(f"{head}_{end}").id for end in ("origin", "distal")] for head in HEADS]
    neutral, _ = measure(model, data, site_ids)
    np.testing.assert_allclose(neutral, data.ten_length, atol=1e-10)
    cases = []
    for peak, active, name in ((0., None, "baseline"), (5., 0, "flex5"),
                               (25., 0, "flex25"), (5., 1, "extend5")):
        rows, result = simulate(model, site_ids, neutral, addresses, 200., peak, active, .005)
        path = output / f"elbow_stop_{name}.csv"
        write_csv(path, ["time_s", *(f"{j}_deg" for j in JOINTS),
                         *(f"{h}_length_change_mm" for h in HEADS),
                         *(f"{h}_drive_N" for h in HEADS),
                         *(f"{h}_TPU_restoring_N" for h in HEADS)], rows)
        cases.append(dict(**result, trace=provenance(path)))
        print(name, "min/max/final elbow",
              *(round(result[f"elbow_{field}_deg"], 5) for field in ("min", "max", "final")),
              "min distance mm", result["minimum_contact_distance_mm"], flush=True)
    save(output / "elbow_local_stop_response.json", dict(
        generator=provenance(Path(__file__).resolve()), local_model=provenance(xml_path),
        force_generator=provenance(ROOT / "tools/simulate_elbow_antagonistic_load.py"),
        stop_reference=provenance(stop_check_path), passive_reference=provenance(passive_path),
        timestep_s=float(model.opt.timestep), trace_interval_s=.005, cases=cases,
        scope="The two existing nominal stop faces are the only contact pair. Negative distance is "
        "numerical penetration. Original provisional masses, frictionless contact law, and assumed "
        "200 N/m linear muscle springs. No full rigid-body, TPU or cable contact certification."))


if __name__ == "__main__":
    main()
