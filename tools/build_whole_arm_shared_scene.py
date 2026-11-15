"""Build the shared arm scene with transferred mounts and original TPU geometry."""

import copy
import struct
import xml.etree.ElementTree as ET

import mujoco
import numpy as np

from audit_s1_series_elbow_integration import find_hits, transforms
from build_channel_first_body_chain import vector
from build_channel_first_inventory import DEST, ROOT, provenance, read
from build_s1_teres_origin_mount import load_component, save
from project_runtime import load_mujoco_model


def main():
    folder = DEST.parent / "whole_arm_first_pass"
    transfer_file = folder / "interface_transfer_audit.json"
    transfer = read(transfer_file)
    if transfer["revised_overlaps_mm3"] or transfer["unassigned"]:
        raise ValueError("The transferred interfaces are still blocked at neutral")
    chain = read(ROOT / transfer["shared_chain"]["file"])
    model_file = ROOT / chain["generated_model"]["file"]
    donor_report = DEST.parent / "integration_review/s1_series_elbow/nominal_spring_spool_seating_k5_coupling.json"
    donor_model = ROOT / read(donor_report)["model"]["file"]
    xml = ET.parse(model_file)
    root = xml.getroot()
    root.set("model", "whole_arm_first_pass")
    asset = root.find("asset")
    bodies = {node.get("name"): node for node in root.findall(".//body")}
    origins = {name: np.asarray(value) * 1000. for name, value in chain["neutral_body_origins_m"].items()}
    origins["thorax"] = np.zeros(3)
    donor = ET.parse(donor_model)
    for parent, name in (("humerus", "FJ1477_carriage"), ("forearm", "FJ1486_distal_pivot"),
                         ("forearm", "FJ1477_distal_pivot")):
        node = copy.deepcopy(donor.find(f".//body[@name='{name}']"))
        for child in list(node):
            if child.tag not in ("joint", "inertial"):
                node.remove(child)
        offset = np.fromstring(node.get("pos", "0 0 0"), sep=" ") * 1000.
        origins[name] = origins[parent] + offset
        bodies[parent].append(node)
        bodies[name] = node

    sources = {**transfer["shared_sources"], **transfer["added_sources"], **transfer["revisions"]}
    owners = transfer["owner_by_source"]
    for row in chain["geom_owners"]:
        key = row["source_key"]
        node = bodies[row["owner"]].find(f"geom[@name='{row['geom']}']")
        bodies[row["owner"]].remove(node)
    for node in list(asset):
        asset.remove(node)

    def add_mesh(key, source, owner, color, group, index):
        path = ROOT / source["file"]
        if provenance(path)["sha256"] != source["sha256"]:
            raise ValueError(f"Invalid mesh source: {key}")
        name = f"whole_arm_mesh_{index}"
        if path.suffix.lower() == ".stl":
            with path.open("rb") as stream:
                stream.seek(80)
                count = struct.unpack("<I", stream.read(4))[0]
            if count > 200000:
                dtype = np.dtype([("normal", "<f4", 3), ("triangle", "<f4", (3, 3)), ("attr", "<u2")])
                triangles = np.frombuffer(path.read_bytes(), dtype=dtype, offset=84)["triangle"]
                vertices, index_map = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
                faces = index_map.reshape(-1, 3) + 1
                preview_assets = folder / "view_assets"
                preview_assets.mkdir(exist_ok=True)
                path = preview_assets / f"{name}.obj"
                with path.open("w", encoding="ascii", newline="\n") as stream:
                    np.savetxt(stream, vertices, fmt="v %.9g %.9g %.9g")
                    np.savetxt(stream, faces, fmt="f %d %d %d")
        ET.SubElement(asset, "mesh", name=name, file=path.relative_to(ROOT).as_posix(), scale=".001 .001 .001")
        offset = np.asarray(source.get("translation_source_mm") or [0., 0., 0.]) - origins[owner]
        ET.SubElement(bodies[owner], "geom", name=name, type="mesh", mesh=name,
                      pos=vector(offset * .001), rgba=color, group=str(group),
                      contype="0", conaffinity="0")

    for index, (key, source) in enumerate(sources.items()):
        color = ".83 .8 .7 1" if key in transfer["shared_sources"] else ".53 .63 .72 1"
        add_mesh(key, source, owners[key], color, 1, index)

    # Existing lattices are shown at neutral. The rigid model remains separate
    # so the meshes cannot be confused with simulated muscle deformation.
    tpu_file = DEST.parent / "integration_review/s1_series_elbow/assembly_FJ1504_local_interface_bidirectional.json"
    tpu = read(tpu_file)["flexible_geometry_sources"]
    preview = ET.Element("body", name="neutral_TPU_reference")
    bodies["thorax"].append(preview)
    bodies["neutral_TPU_reference"] = preview
    origins["neutral_TPU_reference"] = np.zeros(3)
    for index, (key, source) in enumerate(tpu.items(), len(sources)):
        add_mesh(key, source, "neutral_TPU_reference", ".78 .24 .18 .75", 3, index)

    scene = folder / "whole_arm_shared.xml"
    xml.write(scene, encoding="utf-8", xml_declaration=True)
    model = load_mujoco_model(scene)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    for name, origin in origins.items():
        np.testing.assert_allclose(data.body(name).xpos * 1000., origin, atol=1e-7, rtol=0.)
    keys = list(sources)
    solids = [load_component(source) for source in sources.values()]
    pairs = np.triu_indices(len(keys), 1)
    cross_pairs = np.asarray([(i, j) for i, j in zip(*pairs) if owners[keys[i]] != owners[keys[j]]])
    frames = transforms(model, data, {owner: origins[owner] for owner in set(owners.values())})
    neutral_solids = [solid.transform(frames[owners[key]]) for key, solid in zip(keys, solids)]
    hits = find_hits(neutral_solids, keys, cross_pairs)
    if hits:
        raise ValueError(f"Reloaded shared scene is blocked at neutral: {hits}")
    save(folder / "whole_arm_shared_scene_check.json", dict(
        model=provenance(scene), transfer=provenance(transfer_file), donor_kinematics=provenance(donor_model),
        neutral_tpu_reference=provenance(tpu_file), neutral_body_origins_mm={k: v.tolist() for k, v in origins.items()},
        rigid_solid_count=len(sources), visual_tpu_count=len(tpu), geom_count=model.ngeom, nv=model.nv,
        neutral_cross_owner_pairs=len(cross_pairs), neutral_cross_owner_hits=hits,
        scope="Rigid geometry and neutral TPU reference. Group 3 contains fixed neutral TPU previews. "
              "No mesh-contact dynamics, full cable drive, or motion acceptance."))
    print(f"Loaded shared scene: {model.nv} velocity DOFs, {len(sources)} rigid objects, {len(tpu)} TPU references")


if __name__ == "__main__":
    main()
