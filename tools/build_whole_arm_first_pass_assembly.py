"""Assemble the shared K1/K7 spring-relief parts into one review MJCF."""

import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from build_channel_first_inventory import DEST, ROOT, provenance
from build_s1_teres_origin_mount import save
from project_runtime import load_mujoco_model
from study_s1_local_force_geometry import load_reference


def checked_file(source):
    if provenance(ROOT / source["file"])["sha256"] != source["sha256"]:
        raise ValueError(f"Mesh no longer matches its checked source: {source['file']}")
    return source["file"]


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    k1_ref = provenance(folder / "whole_arm_k1_spring_relief_refined.json")
    k7_ref = provenance(folder / "whole_arm_k7_spring_relief.json")
    k1, k7 = map(load_reference, (k1_ref, k7_ref))
    path = load_reference(k1["path"])
    if k7["spring_relief_reference"] != k1_ref or k7["combined"] != path["combined"]:
        raise ValueError("K1 and K7 geometry checks refer to different assemblies")
    combined = load_reference(k7["combined"])
    motion = load_reference(combined["motion"])
    scene = load_reference(motion["scene_check"])
    transfer = load_reference(scene["transfer"])
    support = load_reference(combined["support"])
    sources = load_reference(k1["muscle_meshes"])["tpu_sources"]
    if k7["tpu_sources"]["tpu:deltoid"] != k1["trial"]["tpu:deltoid"]["source"]:
        raise ValueError("K7 did not use the selected deltoid")
    if (any(s["invalid_meshes"] for s in k7["samples"])
            or any(s["spring_contacts_mm3"][key] for s in k7["samples"] for key in ("rigid", "wire", "moving_cuff"))):
        raise ValueError("Spring or TPU geometry checks failed")

    neutral_path = ROOT / path["assembly"]["file"]
    if provenance(neutral_path) != path["assembly"]:
        raise ValueError("Neutral reference changed")
    xml = ET.parse(neutral_path)
    assets = xml.find("asset")
    replacements = {}
    originals = {
        "tpu:deltoid": sources["tpu:deltoid"],
        "elbow:FJ1477_bone_adapter": support["source"],
        "elbow:FJ1477_proximal_rigid_cap": transfer["added_sources"]["elbow:FJ1477_proximal_rigid_cap"],
    }
    for key, old in originals.items():
        new = k1["trial"][key]["source"]
        if not new or new["translation_source_mm"] != old["translation_source_mm"]:
            raise ValueError(f"Part lost its installation frame: {key}")
        meshes = [mesh for mesh in assets.findall("mesh") if mesh.get("file") == checked_file(old)]
        if len(meshes) != 1:
            raise ValueError(f"Cannot resolve the source assembly asset for {key}")
        meshes[0].set("file", checked_file(new))
        replacements[key] = new

    spring = load_reference(k1["spring_definition"])["exports"]["neutral"]["mesh"]
    spring_name = "FJ1477_return_spring_neutral"
    humerus = xml.find(".//body[@name='humerus']")
    position = (np.asarray(spring["translation_source_mm"]) - scene["neutral_body_origins_mm"]["humerus"]) * .001
    ET.SubElement(assets, "mesh", name=spring_name, file=checked_file(spring), scale=".001 .001 .001")
    ET.SubElement(humerus, "geom", name=spring_name, mesh=spring_name, type="mesh", group="2",
                  rgba=".78 .80 .86 1", mass="0", contype="0", conaffinity="0", pos=" ".join(map(str, position)))
    installed = folder / "whole_arm_first_pass_assembled.xml"
    xml.write(installed, encoding="utf-8", xml_declaration=True)
    old_model, model = map(load_mujoco_model, (neutral_path, installed))
    if model.ngeom != old_model.ngeom + 1 or model.nq != old_model.nq or np.count_nonzero(model.geom_group == 3) != 8:
        raise ValueError("Round-trip changed the motion coordinates or lost a mesh")
    for field in ("jnt_type", "jnt_qposadr", "jnt_dofadr", "qpos0"):
        np.testing.assert_array_equal(getattr(model, field), getattr(old_model, field))
    if model.geom(spring_name).bodyid.item() != model.body("humerus").id:
        raise ValueError("The spring is not mounted on the humerus")
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    if not np.isfinite(data.geom_xpos).all():
        raise ValueError("Non-finite assembled geometry")
    save(folder / "whole_arm_first_pass_assembly.json", dict(
        installed=provenance(installed), neutral_reference=path["assembly"], k1=k1_ref, k7=k7_ref,
        replacements=replacements, spring_mesh=spring, guide_points=provenance(folder / "whole_arm_guide_points.json"),
        geometry_count=model.ngeom, qpos_count=model.nq,
        generator=provenance(ROOT / "tools/build_whole_arm_first_pass_assembly.py"),
        scope="Neutral review assembly with a posed return spring, eight TPU meshes and the shared engineered joints. "
              "Exact loaded cuff parts and local wire lines are placed by the K7 replay. The spring is visual here; "
              "deformation, forces and assembly tolerance are separate checks."))
    print(f"Exported and verified {model.ngeom} geometries and {model.nq} position coordinates")


if __name__ == "__main__":
    main()
