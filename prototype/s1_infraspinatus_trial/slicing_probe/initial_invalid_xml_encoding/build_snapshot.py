"""Pack exact source pose plus a Prusa SupportBlocker; never subtract TPU."""
from pathlib import Path
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'mujoco_arm/tools'))
import audit_shoulder_retained_joint_fit as A
import mesh_geom

DEST = ROOT/'prototype/s1_infraspinatus_trial/slicing_probe'
NS = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'


def main():
    manifest = json.loads((DEST/'orientation_manifest.json').read_text(encoding='utf8'))
    trial = json.loads((DEST.parent/'manifest.json').read_text(encoding='utf8'))
    route = trial['active_route_source_mm']
    full = np.array(route['full_wire_path_source_mm'])
    # Extend the blocker over the whole known line, including the actual
    # terminal openings. This is a slicer modifier, not printed geometry.
    guard = A.D.outer_envelope(full, 1.25)
    results = []
    for pose in manifest['orientations']:
        if pose['name'] not in ['paper_z_up', 'paper_z_down']: continue
        source = ROOT/pose['file']
        assert hashlib.sha256(source.read_bytes()).hexdigest() == pose['sha256']
        vertices, faces = mesh_geom.load_obj(str(source))
        affine = np.array(pose['source_to_print_affine'])
        transformed = guard.transform(affine)
        # Prevent the nonprinting volume from enlarging the object bounds
        # used for centering; do not change any printable coordinates.
        lower, upper = vertices.min(axis=0), vertices.max(axis=0)
        box = A.I._manifold().Manifold.cube((upper-lower).tolist()).translate(lower)
        blocker = transformed ^ box
        bv, bf = mesh_geom.manifold_mesh_to_arrays(blocker.to_mesh64())
        model = ET.Element('model', unit='millimeter', attrib={'xmlns': NS, 'xml:lang': 'en-US'})
        ET.SubElement(model, 'metadata', name='slic3rpe:Version3mf').text = '1'
        ET.SubElement(model, 'metadata', name='Application').text = 'PrusaSlicer-2.9.6'
        resources = ET.SubElement(model, 'resources')
        obj = ET.SubElement(resources, 'object', id='1', type='model')
        mesh = ET.SubElement(obj, 'mesh')
        vl, tl = ET.SubElement(mesh, 'vertices'), ET.SubElement(mesh, 'triangles')
        for p in np.vstack([vertices, bv]):
            ET.SubElement(vl, 'vertex', **{axis: format(float(p[i]), '.17g') for i, axis in enumerate('xyz')})
        combined = np.vstack([faces, bf+len(vertices)])
        for face in combined:
            ET.SubElement(tl, 'triangle', **{f'v{i+1}': str(int(face[i])) for i in range(3)})
        ET.SubElement(ET.SubElement(model, 'build'), 'item', objectid='1')
        config = ET.Element('config')
        co = ET.SubElement(config, 'object', id='1', instances_count='1')
        ET.SubElement(co, 'metadata', type='object', key='name', value='ANALYSIS ONLY - full muscle with nonprinting channel blocker')
        for name, volume_type, start, end in [
            ('Current whole TPU muscle', 'ModelPart', 0, len(faces)-1),
            ('Channel keep-open SUPPORT BLOCKER NOT MATERIAL', 'SupportBlocker', len(faces), len(combined)-1)]:
            volume = ET.SubElement(co, 'volume', firstid=str(start), lastid=str(end))
            ET.SubElement(volume, 'metadata', type='volume', key='name', value=name)
            ET.SubElement(volume, 'metadata', type='volume', key='volume_type', value=volume_type)
            ET.SubElement(volume, 'metadata', type='volume', key='matrix', value='1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1')
        target = source.parent/'channel_blocker_ANALYSIS_ONLY.3mf'
        with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED) as package:
            package.writestr('[Content_Types].xml', '<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>')
            package.writestr('_rels/.rels', '<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" Target="/3D/3dmodel.model"/></Relationships>')
            package.writestr('3D/3dmodel.model', ET.tostring(model, encoding='utf8', xml_declaration=True))
            package.writestr('Metadata/Slic3r_PE_model.config', ET.tostring(config, encoding='utf8', xml_declaration=True))
        # Independent readback of the archive's printable triangle block.
        with zipfile.ZipFile(target) as package:
            restored = ET.fromstring(package.read('3D/3dmodel.model'))
            rv = np.array([[float(row.attrib[k]) for k in 'xyz'] for row in restored.findall('.//{'+NS+'}vertex')])
            rf = np.array([[int(row.attrib[k]) for k in ['v1', 'v2', 'v3']] for row in restored.findall('.//{'+NS+'}triangle')])
        assert np.array_equal(rv[:len(vertices)], vertices)
        assert np.array_equal(rf[:len(faces)], faces)
        results.append(dict(orientation=pose['name'], file=target.relative_to(ROOT).as_posix(),
            sha256=hashlib.sha256(target.read_bytes()).hexdigest(), source_pose_sha256=pose['sha256'],
            material_triangles_preserved=len(faces), blocker_triangles=len(bf),
            blocker_radius_mm=1.25, blocker_clipped_to_original_bounds=True,
            volume_types=['ModelPart', 'SupportBlocker'], material_geometry_unchanged=True,
            slicer_roundtrip_type_verified=False))
    report = dict(scope=__doc__, packages=results, format_sources=[
        'https://raw.githubusercontent.com/prusa3d/PrusaSlicer/version_2.9.6/src/libslic3r/Format/3mf.cpp',
        'https://raw.githubusercontent.com/prusa3d/PrusaSlicer/version_2.9.6/src/libslic3r/Model.cpp'],
        physical_print_verified=False, stage_pass=False)
    (DEST/'blocker_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Two source-exact ModelPart + SupportBlocker 3MF probes written.')


if __name__ == '__main__': main()
