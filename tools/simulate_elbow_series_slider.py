"""Drive the local elbow with a free, spring-return FJ1477 carriage."""

from pathlib import Path
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from build_channel_first_body_chain import uniform_geometry_inertia
from simulate_elbow_antagonistic_load import HEADS, JOINTS, write_csv


SLIDE = "FJ1477_series_travel"
SLIDE_STIFFNESS = 73.
SLIDE_DAMPING = .2


def fmt(values):
    return " ".join(str(value) for value in values)


def make_model(parent, guide, check, assembly):
    root = ET.fromstring(parent)
    humerus = root.find(".//body[@name='humerus']")
    source_site = humerus.find("./site[@name='FJ1477_origin']")
    humerus.remove(source_site)
    reference = np.asarray(check["neutral_origins_m"]["humerus"])
    carriage = ET.SubElement(humerus, "body", name="FJ1477_carriage")
    ET.SubElement(carriage, "joint", name=SLIDE, type="slide", axis=fmt(guide["axis_source"]),
                  limited="true", range=fmt([0., guide["travel_mm"] * .001]),
                  solreflimit="-20000000 -6000", solimplimit=".9999 .9999 .0001")
    carriage.append(source_site)
    pieces = []
    source_cap = guide["selected"]["sources"]["moving_cap"]
    fastener_names = [key for key in assembly["sources"] if key.startswith("elbow:FJ1477_proximal_")
                      and key != "elbow:FJ1477_proximal_rigid_cap"]
    for row, density in [(source_cap, 1200.), *((assembly["sources"][key], 7800.) for key in fastener_names)]:
        mass = float(load_component(row).volume()) * 1e-9 * density
        centre, tensor, _ = uniform_geometry_inertia(row, reference, mass)
        pieces.append((mass, centre, tensor))
    mass = sum(piece[0] for piece in pieces)
    centre = sum(m * c for m, c, _ in pieces) / mass
    tensor = np.zeros((3, 3))
    for m, c, inertia in pieces:
        shift = c - centre
        tensor += inertia + m * ((shift @ shift) * np.eye(3) - np.outer(shift, shift))
    ET.SubElement(carriage, "inertial", mass=str(mass), pos=fmt(centre),
                  fullinertia=fmt(tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]))
    actuator = ET.SubElement(root, "actuator")
    for head in HEADS:
        ET.SubElement(actuator, "motor", name=head, tendon=head, gear="-1",
                      ctrllimited="true", ctrlrange="0 100")
    temporary = mujoco.MjModel.from_xml_string(ET.tostring(root, encoding="unicode"))
    data = mujoco.MjData(temporary)
    mujoco.mj_forward(temporary, data)
    for tendon, stiffness, length in zip(root.findall("./tendon/spatial"), (200., 151.53),
                                          data.ten_length, strict=True):
        tendon.set("springlength", str(float(length)))
        tendon.set("stiffness", str(stiffness))
        tendon.set("damping", "1")
    return root, mass, data.ten_length.copy()


def run(model, preload, friction, flex_force, ext_force):
    data = mujoco.MjData(model)
    slide = model.joint(SLIDE)
    model.dof_frictionloss[slide.dofadr.item()] = friction
    # Recompile constants affected by a dof friction change before evaluating.
    mujoco.mj_setConst(model, data)
    addresses = [model.joint(name).qposadr.item() for name in JOINTS]
    slide_pos, slide_dof = slide.qposadr.item(), slide.dofadr.item()
    timestep = model.opt.timestep
    stroke = model.jnt_range[slide.id, 1]
    rest = model.tendon_lengthspring[:, 0].copy()
    rows = []
    upper_peak = lower_peak = 0.
    upper_steps = lower_steps = 0
    stop_min = 0.
    for index in range(round(4. / timestep) + 1):
        mujoco.mj_step1(model, data)
        time = index * timestep
        pos, vel = float(data.qpos[slide_pos]), float(data.qvel[slide_dof])
        restoring = max(0., preload + SLIDE_STIFFNESS * pos + SLIDE_DAMPING * vel)
        data.qfrc_applied[slide_dof] = -restoring
        flex = flex_force * np.interp(time, [0., .5, 1.5, 2.5, 3., 4.], [0., 0., 1., 1., 0., 0.])
        ext = ext_force * np.interp(time, [0., 2.5, 2.75, 3.2, 3.5, 4.], [0., 0., 1., 1., 0., 0.])
        data.ctrl[:] = [flex, ext]
        upper_steps += int(pos > stroke)
        lower_steps += int(pos < 0.)
        upper_peak = max(upper_peak, pos - stroke)
        lower_peak = max(lower_peak, -pos)
        if data.ncon:
            stop_min = min(stop_min, float(min(data.contact.dist[:data.ncon])))
        if index % round(.0025 / timestep) == 0:
            rows.append([time, *np.rad2deg(data.qpos[addresses]), pos * 1000.,
                         *((data.ten_length - rest) * 1000.), flex, ext, restoring])
        if time < 4. - timestep / 2.:
            mujoco.mj_step2(model, data)
    values = np.asarray(rows)
    if not np.isfinite(values).all() or any(w.number for w in data.warning):
        raise RuntimeError("Invalid free-slider response")
    return rows, dict(
        slider_preload_N=preload, rail_frictionloss_N=friction, flexor_peak_force_N=flex_force,
        extensor_peak_force_N=ext_force, elbow_peak_deg=float(values[:, 1].max()),
        elbow_final_deg=float(values[-1, 1]), carriage_peak_mm=float(values[:, 4].max()),
        carriage_at_loading_start_mm=float(values[200, 4]), carriage_final_mm=float(values[-1, 4]),
        peak_FJ1477_span_growth_mm=float(values[:, 6].max()),
        peak_series_return_force_N=float(values[:, 9].max()),
        upper_limit_overshoot_mm=upper_peak * 1000., lower_limit_overshoot_mm=lower_peak * 1000.,
        upper_limit_overshoot_steps=upper_steps, lower_limit_overshoot_steps=lower_steps,
        elbow_stop_contact_min_distance_mm=stop_min * 1000.,
        max_saved_angle_step_deg=float(np.max(abs(np.diff(values[:, 1:4], axis=0)))),
        max_saved_carriage_step_mm=float(np.max(abs(np.diff(values[:, 4])))),
    )


def main():
    output = DEST.parent / "elbow_local_muscle_load"
    guide_file = output / "elbow_series_guide.json"
    audit_file = output / "elbow_series_guide_audit.json"
    stop_file = output / "elbow_local_stop_response.json"
    guide, audit, stop = map(read, (guide_file, audit_file, stop_file))
    stop_model = ROOT / stop["local_model"]["file"]
    assembly_file = ROOT / guide["stop_assembly"]["file"]
    if (audit["design"] != provenance(guide_file) or audit["first_cycle_failure"] is not None
            or stop["local_model"] != provenance(stop_model)
            or guide["stop_assembly"] != provenance(assembly_file)):
        raise RuntimeError("Series guide or stop model differs from the audited geometry")
    passive_file = DEST / "elbow_passive_model/model_check.json"
    root, mass, reference_span = make_model(stop_model.read_text(encoding="utf-8"), guide,
                                            read(passive_file), read(assembly_file))
    xml = ET.tostring(root, encoding="unicode")
    model_file = output / "elbow_series_free_slider.xml"
    model_file.write_text(xml, encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(xml)
    if model.nv != 4 or model.ntendon != 2:
        raise RuntimeError("Expected three original hinges and one added slide")
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    np.testing.assert_allclose(data.ten_length, reference_span, atol=1e-10)
    np.testing.assert_allclose(data.site("FJ1477_origin").xpos * 1000.,
                               next(port["source_position_mm"] for port in read(passive_file)["ports"]
                                    if port["head"] == "FJ1477" and port["role"] == "origin"), atol=1e-6)
    cases = []
    settings = ((0., 0., 25., 0.), (.25, 0., 25., 0.), (.25, .15, 25., 0.),
                (.25, .5, 25., 0.), (.25, .15, 25., 5.))
    for index, params in enumerate(settings):
        rows, metrics = run(model, *params)
        trace_path = output / f"elbow_series_free_slider_case{index}.csv"
        write_csv(trace_path, ["time_s", *(f"{joint}_deg" for joint in JOINTS),
                              "carriage_mm", *(f"{head}_span_change_mm" for head in HEADS),
                              "flexor_force_N", "extensor_force_N", "series_return_force_N"], rows)
        cases.append(dict(**metrics, trace=provenance(trace_path)))
        print("case", index, "elbow peak/final", metrics["elbow_peak_deg"], metrics["elbow_final_deg"],
              "carriage peak/final", metrics["carriage_peak_mm"], metrics["carriage_final_mm"], flush=True)
    save(output / "elbow_series_free_slider.json", dict(
        generator=provenance(Path(__file__).resolve()), local_model=provenance(model_file),
        guide=provenance(guide_file), guide_audit=provenance(audit_file), stop=provenance(stop_file),
        carriage_mass_kg=mass, return_stiffness_N_per_m=SLIDE_STIFFNESS,
        return_damping_Ns_per_m=SLIDE_DAMPING, muscle_stiffness_N_per_m=[200., 151.53],
        muscle_damping_Ns_per_m=1., cases=cases,
        scope="Unilateral imposed force along endpoint chords; active steel routing and flexibility are "
        "unresolved. The slider and three hinges are free to move. Slider travel bounds use software "
        "limits, and the original stop face pair is retained at the elbow. Passive springs, damping "
        "and rail frictionloss are hypotheses. Return spring is applied by the simulation script. "
        "No collision solver for added guide meshes, physical spring design, or TPU validation."))


if __name__ == "__main__":
    main()
