r"""Generate first-pass PSM geometry inside Rhino from the CAD manifest.

Run this file from Rhino 8 with the command:
    _-RunPythonScript "D:\XJTLU\PSP\mujoco_arm\grasshopper_interface\generate_psm_rhino.py"

The script is intentionally compatible with Rhino's IronPython 2.7 runtime.
"""

import io
import json
import math
import os
import traceback

import Rhino
import rhinoscriptsyntax as rs


ROOT = r"D:\XJTLU\PSP\mujoco_arm"
MANIFEST = os.path.join(ROOT, "exports", "grasshopper_muscles_cad.json")
OUTPUT_DIR = os.path.join(ROOT, "exports", "rhino_psm")
LOG_PATH = os.path.join(OUTPUT_DIR, "generate_psm.log")


def log(message):
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    stream = io.open(LOG_PATH, "a", encoding="utf-8")
    try:
        stream.write(unicode(message) + u"\n")
    finally:
        stream.close()


def as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def add_pipe(curve, radius, layer):
    objects = as_list(
        rs.AddPipe(curve, [0.0, 1.0], [radius, radius], 0, 2)
    )
    if not objects:
        raise RuntimeError("Rhino failed to create a pipe.")
    rs.ObjectLayer(objects, layer)
    return objects


def add_lattice(base_plane, length, width, cell_size, branch_radius, layer):
    longitudinal = max(3, int(math.ceil(length / cell_size)) + 1)
    transverse = 3
    segment = length / float(longitudinal - 1)
    transverse_spacing = width / float(transverse - 1)
    v = base_plane.YAxis
    w = base_plane.ZAxis
    origin = base_plane.Origin

    def node(i, j, sign):
        return (
            origin
            + base_plane.XAxis * (i * segment)
            + v * ((j - (transverse - 1) / 2.0) * transverse_spacing)
            + w * (sign * width * 0.5)
        )

    for i in range(longitudinal - 1):
        sign = 1 if i % 2 == 0 else -1
        for j in range(transverse):
            curve = rs.AddLine(node(i, j, sign), node(i + 1, j, -sign))
            if curve is not None:
                add_pipe(curve, branch_radius, layer)

    for i in range(longitudinal):
        sign = 1 if i % 2 == 0 else -1
        for j in range(transverse - 1):
            curve = rs.AddLine(node(i, j, sign), node(i, j + 1, sign))
            if curve is not None:
                add_pipe(curve, branch_radius * 0.85, layer)


def main():
    if os.path.exists(LOG_PATH):
        os.remove(LOG_PATH)
    log("manifest: " + MANIFEST)
    if not os.path.exists(MANIFEST):
        raise RuntimeError("Manifest not found: " + MANIFEST)

    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    stream = io.open(MANIFEST, "r", encoding="utf-8")
    try:
        muscles = json.loads(stream.read())
    finally:
        stream.close()

    rs.EnableRedraw(False)
    try:
        all_objects = []
        for muscle in muscles:
            muscle_id = muscle["id"]
            layer = "PSM::" + muscle_id
            if not rs.IsLayer(layer):
                rs.AddLayer(layer)

            origin = rs.CreatePoint(*muscle["origin_world"])
            insertion = rs.CreatePoint(*muscle["insertion_world"])
            path_points = [origin]
            for via in muscle.get("via_world", []):
                path_points.append(rs.CreatePoint(*via["position"]))
            path_points.append(insertion)
            curve = rs.AddPolyline(path_points)
            cable_layer = layer + "::Cable"
            rs.ObjectLayer(curve, cable_layer)
            all_objects.extend(add_pipe(curve, 0.00075, cable_layer))

            axis = insertion - origin
            length = axis.Length
            if length < 1e-6:
                continue
            axis.Unitize()
            reference = Rhino.Geometry.Vector3d(0, 0, 1)
            if abs(axis * reference) > 0.95:
                reference = Rhino.Geometry.Vector3d(0, 1, 0)
            side = Rhino.Geometry.Vector3d.CrossProduct(axis, reference)
            side.Unitize()
            normal = Rhino.Geometry.Vector3d.CrossProduct(side, axis)
            normal.Unitize()
            plane = Rhino.Geometry.Plane(origin, axis, side)

            lattice = muscle["psm_parameters"]
            x_m = lattice["x_mm"] * 0.001
            y_m = lattice["y_mm"] * 0.001
            p_m = lattice["p_mm"] * 0.001
            lattice_layer = layer + "::Lattice"
            add_lattice(
                plane,
                length,
                y_m * 2.0,
                x_m,
                p_m * 0.5,
                lattice_layer,
            )
            all_objects.extend(as_list(rs.ObjectsByLayer(lattice_layer)))

            end_layer = layer + "::Ends"
            for point in (origin, insertion):
                sphere = rs.AddSphere(
                    point, max(y_m * 0.75, p_m * 1.5)
                )
                if sphere is not None:
                    rs.ObjectLayer(sphere, end_layer)
                    all_objects.append(sphere)

        if not all_objects:
            raise RuntimeError("No PSM geometry was created.")

        rs.EnableRedraw(True)
        rs.ZoomExtents()
        selected = rs.SelectObjects(all_objects)
        if not selected:
            raise RuntimeError("Unable to select generated PSM geometry.")

        model_path = os.path.join(OUTPUT_DIR, "right_arm_psm.3dm")
        stl_path = os.path.join(OUTPUT_DIR, "right_arm_psm.stl")
        rs.Command('_-SaveAs "' + model_path + '" _Enter', False)
        rs.Command('_-Export "' + stl_path + '" _Enter _Enter', False)
        log("Rhino model: " + model_path)
        log("STL export: " + stl_path)
        if os.environ.get("PSM_RHINO_EXIT") == "1":
            rs.Command("_Exit", False)
    finally:
        rs.EnableRedraw(True)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        log(traceback.format_exc())
        raise
