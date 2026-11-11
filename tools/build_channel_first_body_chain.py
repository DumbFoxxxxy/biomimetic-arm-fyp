"""Build an explicit ideal joint-chain reference for the channel-first model."""
import argparse
import sys
import xml.etree.ElementTree as ET

import numpy as np

from build_s1_teres_origin_mount import ROOT, mesh_geom, read, save
from build_channel_first_inventory import provenance


def vector(values):
    return ' '.join(format(float(value), '.12g') for value in values)


def uniform_geometry_inertia(source, origin, mass):
    path = ROOT/source['file']
    vertices, faces = (mesh_geom.load_obj(str(path)) if path.suffix.lower() == '.obj'
                       else mesh_geom.read_stl(str(path)))
    vertices = (vertices+np.array(source.get('translation_source_mm') or [0., 0., 0.]))*.001-origin
    triangle = vertices[faces]
    a, b, c = triangle[:, 0], triangle[:, 1], triangle[:, 2]
    volumes = np.einsum('ij,ij->i', a, np.cross(b, c))/6.
    volume = volumes.sum()
    assert volume > 0.
    sums = triangle.sum(axis=1)
    centre = np.sum(volumes[:, None]*sums/4., axis=0)/volume
    second = np.einsum('n,nij->ij', volumes,
                        np.einsum('nki,nkj->nij', triangle, triangle)
                        +np.einsum('ni,nj->nij', sums, sums))/20./volume
    covariance = second-np.outer(centre, centre)
    inertia = mass*(np.trace(covariance)*np.eye(3)-covariance)
    assert np.linalg.eigvalsh(inertia).min() > 0.
    return centre, inertia, float(volume)


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry-mass', action='store_true')
    parser.add_argument('--unified', action='store_true', help='Use selected unified skeleton in its own output folder')
    parser.add_argument('--reference', help='Explicit alternative unified assembly; requires --unified')
    args = parser.parse_args()

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    reference_file = folder/'thorax_arm_reference.json'
    if args.unified:
        pointer = read(ROOT/'design/CHANNEL_FIRST_SKELETON.json')
        reference_file = ROOT/(args.reference or pointer['assembly']['file'])
        if not args.reference:
            assert provenance(reference_file) == pointer['assembly']
        folder = reference_file.parent
    reference = read(reference_file)
    definition_file = ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/definition.json'
    definition = read(definition_file)
    motion_file = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus/anatomical_elbow_rigid_motion.json'
    elbow = read(motion_file)
    centres = {key: np.asarray(value)*.001 for key, value in definition['centres_source_mm'].items()}
    assert np.allclose(centres['elbow']*1000., elbow['centre_source_mm'])
    sources = dict(reference['sources'])
    joint_owners = {}
    if args.unified:
        sources.update(reference['joint_context_sources'])
        joint_owners = {key: ('thorax' if key.startswith('SC:') else
                            'forearm' if key == 'elbow:elbow_bearing_bush' else
                            'humerus' if key.startswith('elbow:') else 'scapula')
                        for key in reference['joint_context_sources']}
    # Mass and scalar inertia are provisional identification parameters.
    definitions = [
        ('clavicle', 'thorax', 'SC', 'ball', .20, .0005),
        ('scapula', 'clavicle', 'AC', 'ball', .30, .0010),
        ('humerus', 'scapula', 'GH', 'ball', .60, .0060),
        ('forearm', 'humerus', 'elbow', 'hinge', .40, .0040)]
    root = ET.Element('mujoco', model='channel_first_ideal_body_chain')
    ET.SubElement(root, 'compiler', angle='radian')
    ET.SubElement(root, 'option', timestep='.001', gravity='0 0 -9.81', integrator='implicitfast')
    asset = ET.SubElement(root, 'asset')
    world = ET.SubElement(root, 'worldbody')
    thorax = ET.SubElement(world, 'body', name='thorax')
    body_nodes = {'thorax': thorax}
    body_origins = {'thorax': np.zeros(3)}
    body_report = []
    for name, parent, joint, kind, mass, inertia in definitions:
        origin = centres[joint]
        node = ET.SubElement(body_nodes[parent], 'body', name=name,
                             pos=vector(origin-body_origins[parent]))
        attributes = dict(name=joint, type=kind, damping='.05')
        if kind == 'hinge':
            attributes.update(axis=vector(elbow['axis_source']), limited='true', range=vector(np.radians([-5, 150])))
        ET.SubElement(node, 'joint', **attributes)
        if args.geometry_mass:
            com, tensor, volume = uniform_geometry_inertia(sources[f'bone:{name}'], origin, mass)
            full = tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
            ET.SubElement(node, 'inertial', pos=vector(com), mass=str(mass), fullinertia=vector(full))
        else:
            com, tensor, volume = np.zeros(3), np.eye(3)*inertia, None
            ET.SubElement(node, 'inertial', pos='0 0 0', mass=str(mass), diaginertia=vector([inertia]*3))
        body_nodes[name], body_origins[name] = node, origin
        body_report.append(dict(body=name, parent=parent, joint=joint, type=kind,
                                joint_centre_source_mm=(origin*1000.).tolist(),
                                assumed_mass_kg=mass, centre_of_mass_local_m=com.tolist(),
                                inertia_tensor_kg_m2=tensor.tolist(), geometry_volume_m3=volume))
    geom_rows = []
    for index, (key, source) in enumerate(sources.items()):
        path = ROOT/source['file']
        assert provenance(path)['sha256'] == source['sha256']
        owner = joint_owners.get(key, key.removeprefix('bone:') if key.startswith('bone:') else 'thorax')
        assert owner in body_nodes
        name = f'mesh_{index:03d}'
        ET.SubElement(asset, 'mesh', name=name, file=str(path), scale='.001 .001 .001')
        translation = np.asarray(source.get('translation_source_mm') or [0., 0., 0.])*.001
        position = translation-body_origins[owner]
        ET.SubElement(body_nodes[owner], 'geom', name=name, type='mesh', mesh=name, pos=vector(position),
                      contype='0', conaffinity='0', group='1',
                      rgba='.88 .84 .74 1' if owner != 'thorax' else '.65 .7 .75 1')
        geom_rows.append(dict(geom=name, source_key=key, owner=owner))
    stem = 'geometry_mass_body_chain' if args.geometry_mass else 'ideal_body_chain'
    xml_file = folder/f'{stem}.xml'
    ET.ElementTree(root).write(xml_file, encoding='utf-8', xml_declaration=True)
    model = load_mujoco_model(xml_file)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    positions = {}
    for row in body_report:
        body_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, row['body'])
        expected = np.array(row['joint_centre_source_mm'])*.001
        assert np.allclose(data.xpos[body_id], expected, atol=1e-10)
        positions[row['body']] = data.xpos[body_id].tolist()
    initial = data.qpos.copy()
    for _ in range(100):
        mujoco.mj_step(model, data)
    assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
    save(folder/f'{stem}_check.json', dict(
        reference=provenance(reference_file), anatomical_definition=provenance(definition_file),
        elbow_reference=provenance(motion_file), generated_model=provenance(xml_file),
        bodies=body_report, geom_owners=geom_rows, neutral_body_origins_m=positions,
        nq=model.nq, nv=model.nv, geom_count=model.ngeom,
        gravity_smoke_test_seconds=.1, initial_qpos=initial.tolist(), final_qpos=data.qpos.tolist(),
        final_qvel=data.qvel.tolist(), finite_gravity_steps=True,
        limitations=[
            'SC and AC are ideal spherical research joints, replacing historical fixed bench conditions',
            ('geometry uses the selected unified bone candidate and explicitly assigned GH/elbow context'
             if args.unified else 'geometry uses the explicit thorax_arm_reference arm sources, not a merged claim of all S1 revisions'),
            'thorax is mathematically grounded; physical base and thoracic assembly are pending',
            ('mass and damping are assumed; centre of mass and inertia use uniform-density bone geometry'
             if args.geometry_mass else 'mass, inertia and damping are provisional assumed parameters'),
            'geoms are visual only; anatomical mesh convex hulls are not collision representations',
            'no muscles, cables, guide reactions, or range acceptance yet']))
    print('ideal chain:', model.nv, 'velocity DOFs,', model.ngeom, 'visual geoms; finite gravity steps')


if __name__ == '__main__':
    main()
