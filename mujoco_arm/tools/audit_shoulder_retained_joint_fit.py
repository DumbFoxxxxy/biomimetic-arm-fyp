"""Compare actual retained GH hardware with the existing compact TPU candidate.

This is a zero-pose compatibility experiment, not an adopted assembly. The G1
meshes include earlier AC work; clashes involving that context are retained.
All OBJ/STL inputs are hashed, checked for closed topology and positive volume,
and explicitly transformed into source millimetres before intersections.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frames
import interference as I
import mesh_geom
import design_shoulder_port_paths as D

ROOT = frames.PROJECT_ROOT
DESIGN = frames.OUTPUT_DIR.parent / 'shoulder_design'
G1 = (frames.OUTPUT_DIR / 'joint_model/SC_connection/AC_compact/GH_connection')
CURRENT = DESIGN / 's1_port_path_revision/whole_shoulder_compact_distal'
DEST = DESIGN / 's1_retained_joint_comparison'


def hashed(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_checked(path, translation=None, expected=None):
    if expected is not None and hashed(path) != expected:
        raise ValueError('Changed input: ' + str(path))
    vertices, faces = (mesh_geom.load_obj(str(path)) if path.suffix == '.obj'
                       else mesh_geom.read_stl(str(path)))
    vertices, faces = mesh_geom.weld(vertices, faces)
    _, boundary, nonmanifold = mesh_geom.edge_census(faces)
    if boundary or nonmanifold:
        raise ValueError(f'Open/nonmanifold input {path}: {boundary}/{nonmanifold}')
    if translation is not None:
        vertices = vertices + np.asarray(translation)
    solid = I._manifold().Manifold(
        I._manifold().Mesh64(vertices, faces.astype(np.uint32)))
    if solid.volume() <= 0 or str(solid.status()).split('.')[-1] != 'NoError':
        raise ValueError('Invalid or nonpositive solid: ' + str(path))
    metadata = dict(file=path.relative_to(ROOT).as_posix(), sha256=hashed(path),
                    translation_source_mm=None if translation is None else list(translation),
                    boundary_edges=boundary, nonmanifold_edges=nonmanifold,
                    volume_mm3=float(solid.volume()))
    return solid, vertices[faces], metadata


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    assembly_file = CURRENT / 'actual_assembly_check.json'
    routes_file = DESIGN / 's1_port_path_revision/routes_compact_distal_candidate.json'
    joint_file = G1 / 'assembly.json'
    assembly = json.loads(assembly_file.read_text(encoding='utf8'))
    routes = json.loads(routes_file.read_text(encoding='utf8'))['routes']
    joint = json.loads(joint_file.read_text(encoding='utf8'))
    delta = np.asarray(joint['centre_source_mm']) - frames.joint_centres('source')['GH']
    if np.linalg.norm(delta) > 1e-6:
        raise ValueError('Joint centres differ; explicit registration is required')
    retained, current, source, tris = {}, {}, {}, {}
    invalid_retained = []
    for key in joint['components']:
        path = G1 / (key + '.obj')
        try:
            solid, tri, meta = load_checked(path, frames.SHIFT)
        except ValueError as exc:
            invalid_retained.append(dict(part=key, file=path.relative_to(ROOT).as_posix(),
                sha256=hashed(path), reason=str(exc)))
            print('INVALID retained input', key, str(exc), flush=True)
            continue
        retained[key] = solid
        tris[key] = tri
        source['retained:' + key] = meta
    for key, row in assembly['sources'].items():
        solid, tri, meta = load_checked(ROOT / row['file'],
                                       row['translation_source_mm'], row['sha256'])
        current[key] = solid
        tris[key] = tri
        source[key] = meta
    replacements = {'bone:humerus': 'humerus_GH_ball',
                    'bone:scapula': 'scapula_AC_GH_sockets'}
    pairs = []
    for name, solid in retained.items():
        hits = []
        for key, old in current.items():
            # The two alternate bone solids replace rather than coexist with
            # the anatomical ones. No other bone or attachment pair is skipped.
            if key in replacements:
                continue
            volume = I.intersection_volume(solid, old)
            if volume > 1e-9:
                hits.append(dict(current=key, volume_mm3=volume,
                                 review='compliant contact requires deformation review'
                                 if key.startswith('tpu:') else 'rigid conflict'))
        pairs.extend(dict(retained=name, **hit) for hit in hits)
        print(name, 'cross-candidate overlaps',
              [(h['current'], round(h['volume_mm3'], 5)) for h in hits], flush=True)
    # Separate inherited AC context from material added for the GH socket.
    ac_path = G1.parent / 'scapula_AC_socket.obj'
    ac, _, ac_meta = load_checked(ac_path, frames.SHIFT)
    source['context:scapula_AC_socket'] = ac_meta
    gh_added = retained['scapula_AC_GH_sockets'] - ac
    localization = []
    for key, candidate in current.items():
        if not key.startswith(('tpu:', 'cap:')):
            continue
        intersection = gh_added ^ candidate
        volume = float(intersection.volume())
        if volume > 1e-9:
            localization.append(dict(current=key, GH_added_overlap_mm3=volume,
                source_bounds_mm=list(intersection.bounding_box())))
    seats = json.loads((DESIGN / 's1_end_pads.json').read_text(encoding='utf8'))['seats']
    owners = {(label.split('/')[0], s['end']): s['nearest_bone']
              for s in seats for label in s['heads']}
    attachment_rows = []
    for route in routes:
        for end, port in route['ports'].items():
            bone = owners[route['head'], end]
            original = 'bone:' + bone
            replacement = replacements.get(original)
            target_tri = tris.get(replacement) if replacement else tris[original]
            seat = np.asarray(port['bone_seat_source_mm'])
            old_dist = float(mesh_geom.point_triangle_distance(seat, tris[original]).min())
            distance = (float(mesh_geom.point_triangle_distance(seat, target_tri).min())
                        if target_tri is not None else None)
            attachment_rows.append(dict(head=route['head'], end=end, bone=bone,
                seat_source_mm=seat.tolist(), original_surface_distance_mm=old_dist,
                retained_surface_distance_mm=distance,
                distance_change_mm=None if distance is None else distance-old_dist))
    wires = []
    for row in routes:
        wire = D.outer_envelope(np.asarray(row['full_wire_path_source_mm']), .3)
        hits = [dict(retained=k, volume_mm3=I.intersection_volume(wire, s))
                for k, s in retained.items()]
        wires.append(dict(head=row['head'], nominal_wire_radius_mm=.3,
                          intersections=hits,
                          hits_over_0_001_mm3=[h for h in hits if h['volume_mm3'] > .001]))
        print('wire', row['head'], wires[-1]['hits_over_0_001_mm3'], flush=True)
    report = dict(scope='zero-pose cross-candidate compatibility only; NOT an adopted S1',
        input_sha256={p.relative_to(ROOT).as_posix(): hashed(p)
                      for p in [assembly_file, routes_file, joint_file,
                                DESIGN/'s1_end_pads.json', Path(__file__).resolve()]},
        source_frame='BodyParts3D source millimetres; G1 display + frames.SHIFT',
        centres_difference_mm=delta.tolist(), sources=source,
        replaced_bones=replacements, invalid_retained_inputs=invalid_retained,
        checked_retained_parts=len(retained), expected_retained_parts=len(joint['components']),
        cross_candidate_intersections=pairs, GH_added_conflict_localization=localization,
        attachment_surface_distances=attachment_rows, wire_probes=wires,
        not_verified=['joint adoption or new CAD', 'moving and combined paths',
            'TPU deformation and friction', 'fastening, retention strength and assembly access',
            'G1 includes old AC socket modifications, not an isolated new GH insert'],
        stage_pass=False, actual_drive_range_deg=None)
    (output/'compatibility.json').write_text(json.dumps(report, indent=2)+'\n',
                                           encoding='utf8', newline='\n')
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=DEST)
    run(parser.parse_args().output_dir)
