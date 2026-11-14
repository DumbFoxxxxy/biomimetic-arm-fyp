"""Measure local tendon excursions and cross-joint moment arms on a checked trace."""

import mujoco
import numpy as np
import xml.etree.ElementTree as ET

from audit_s1_series_elbow_timed_muscle_geometry import checked_trace
from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from simulate_elbow_antagonistic_load import write_csv
from study_s1_local_force_geometry import load_reference


def main():
    folder = DEST.parent / "integration_review/s1_series_elbow"
    stem = "nominal_spring_spool_seating_k5"
    flexible_ref = provenance(folder / f"{stem}_flexible_audit.json")
    flexible = load_reference(flexible_ref)
    rigid = load_reference(flexible["spool_rigid_audit"])
    trace = checked_trace(flexible["trace"])
    if (rigid["trace"] != flexible["trace"] or rigid["findings"] or flexible["findings"]
            or len(trace) != rigid["sampled_frames"] or len(trace) != flexible["nominal_rigid_and_wire_frames"]):
        raise ValueError("Coupling study requires matching, clear geometry reports")
    model_path = ROOT / flexible["model"]["file"]
    if provenance(model_path) != flexible["model"]:
        raise ValueError("Drive model changed")
    xml = model_path.read_text(encoding="utf-8")
    model = mujoco.MjModel.from_xml_string(xml)
    paths = {row.get("name"): [point.get("site") for point in row]
             for row in ET.fromstring(xml).find("tendon")}
    data = mujoco.MjData(model)
    mujoco.mj_kinematics(model, data)
    anchors = {name: [dict(site=site, body=model.body(model.site_bodyid[model.site(site).id]).name,
                           reference_xyz_mm=(1000. * data.site(site).xpos).tolist())
                      for site in path] for name, path in paths.items()}
    gh = model.joint("GH").dofadr.item()
    elbow = model.joint("elbow").dofadr.item()
    slider = model.joint("FJ1477_series_travel").dofadr.item()
    dofs = (gh, gh + 1, gh + 2, elbow, slider)
    names = ("gh_x_mm", "gh_y_mm", "gh_z_mm", "elbow_mm", "slider_mm_per_mm")
    epsilon = 1e-6
    rows = []
    samples = []
    for frame in trace:
        qpos = np.array([float(frame[f"qpos_{i}"]) for i in range(model.nq)])
        data.qpos[:] = qpos
        mujoco.mj_forward(model, data)
        lengths = data.ten_length.copy() * 1000.
        gradient = []
        for dof in dofs:
            direction = np.zeros(model.nv)
            direction[dof] = 1.
            for sign in (-1., 1.):
                data.qpos[:] = qpos
                mujoco.mj_integratePos(model, data.qpos, direction, sign * epsilon)
                mujoco.mj_forward(model, data)
                if sign < 0.:
                    minus = data.ten_length.copy()
            gradient.append((minus - data.ten_length) / (2. * epsilon) * (1. if dof == slider else 1000.))
        values = np.column_stack((lengths, np.array(gradient).T))
        samples.append(values)
        for tendon, row in enumerate(values):
            rows.append([frame["time_s"], frame["elbow_deg"], model.tendon(tendon).name, *row])
    values = np.array(samples)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite local tendon measurements")
    channels = {}
    for tendon in range(model.ntendon):
        series = values[:, tendon]
        channels[model.tendon(tendon).name] = dict(
            equivalent_sites=anchors[model.tendon(tendon).name],
            excursion_mm=float(np.ptp(series[:, 0])),
            end_delta_mm=float(series[-1, 0] - series[0, 0]),
            gh_moment_arm_norm_peak_mm=float(np.max(np.linalg.norm(series[:, 1:4], axis=1))),
            derivatives={key: [float(series[:, col].min()), float(series[:, col].max())]
                         for col, key in enumerate(names, start=1)},
        )
    output = folder / f"{stem}_coupling.csv"
    write_csv(output, ["time_s", "elbow_deg", "tendon", "equivalent_length_mm", *names], rows)
    save(folder / f"{stem}_coupling.json", dict(
        generator=provenance(ROOT / "tools/audit_s1_series_elbow_coupling.py"),
        flexible_audit=flexible_ref, trace=flexible["trace"], model=flexible["model"],
        sampled_frames=len(trace), channels=channels, table=provenance(output),
        scope="Negative partial derivatives of equivalent tendon length with all other coordinates fixed. "
              "GH directions are MuJoCo ball-joint tangent coordinates, angular derivatives in mm/rad. "
              "Slider derivative is dimensionless. Model has a fixed shoulder girdle. "
              "Anatomical cable wrap, actuator-side routing, contact and structural cross-talk are excluded.",
    ))
    print(f"Measured {len(trace)} frames for {model.ntendon} tendons")


if __name__ == "__main__":
    main()
