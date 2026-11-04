"""Bench-only integral proximal saddle; does not revise the shoulder assembly.

Extends the actual front support backwards to a vice block. Retains the
two M3 cap axes, provides rear nut/tool wells, and leaves the distal end free.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'mujoco_arm/tools'), str(ROOT)]
import audit_shoulder_retained_joint_fit as A
import mesh_geom
import shoulder_pads

DEST = ROOT/'prototype/s1_infraspinatus_trial/bench_adapter_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')


def export(solid, path, origin):
    local = solid.translate(-origin)
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(local.to_mesh64())
    mesh_geom.write_stl(str(path), vertices, faces, header='Bench-only saddle; no strength or physical drive verification')
    actual, _, metadata = A.load_checked(path, origin.tolist(), sha(path))
    components = [float(piece.volume()) for piece in actual.decompose()]
    if len(components) != 1 or components[0] <= 0:
        raise ValueError('Adapter must remain a single positive solid, without automatic repairs')
    metadata.update(components_mm3=components,
                    roundtrip_relative_volume_error=abs(actual.volume()-solid.volume())/solid.volume())
    if metadata['roundtrip_relative_volume_error'] > 1e-5:
        raise ValueError('Unacceptable export roundtrip')
    return actual, metadata


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    report_file = ROOT/selection['assembly_report']
    if sha(report_file) != selection['assembly_sha256']:
        raise ValueError('Stale current assembly')
    parent = read(report_file)
    trial_file = DEST.parent/'assembly_subset.json'
    trial = read(trial_file)
    route_file = ROOT/read(DEST.parent/'manifest.json')['active_route_source']
    route = next(r for r in read(route_file)['routes'] if r['head'] == 'FJ1500')
    port = route['ports']['origin']
    origin = np.array(port['inner_face_source_mm'])
    basis = shoulder_pads.pad_basis(np.array(port['normal_source']))
    normal = basis[:, 2]
    solids, vertices = {}, {}
    for key, source in trial['sources'].items():
        solid, triangles, _ = A.load_checked(ROOT/source['file'], source['translation_source_mm'], source['sha256'])
        if source['sha256'] != parent['sources'][key]['sha256']:
            raise ValueError('Trial differs from current candidate')
        solids[key], vertices[key] = solid, triangles.reshape(-1, 3)
    mf = A.I._manifold()
    # All construction coordinates are in the real proximal cap basis.
    stem = mf.Manifold.cylinder(37.7, 11., 11., 96).translate([0., 0., -42.])
    grip = mf.Manifold.cube([34., 28., 18.]).translate([-17., -14., -56.])
    addition = (stem+grip).transform(np.c_[basis, origin])
    # Relieve only the added bracket, never the original clamp surface.
    tpu = solids['tpu:infraspinatus']
    for delta in [np.zeros(3)]+[basis[:, i]*sign*.3 for i in range(3) for sign in [-1, 1]]:
        addition = addition-tpu.translate(delta)
    wire_guard = A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), 1.2)
    addition = addition-wire_guard
    wells = []
    for index, sign in enumerate([-1, 1]):
        washer_key = f'origin:FJ1500/back_washer_{index}'
        local_washer = (vertices[washer_key]-origin)@basis
        bearing_z = float(local_washer[:, 2].max())+.02
        # The existing support already contains its own 3.4mm holes.
        # A 3.5mm bore in the ADDED material avoids re-cutting a coincident
        # faceted inherited bore and keeps the existing interface untouched.
        bore = mf.Manifold.cylinder(70., 1.75, 1.75, 64).translate([sign*7., 0., -60.])
        access = mf.Manifold.cylinder(bearing_z+60., 4.5, 4.5, 96).translate([sign*7., 0., -60.])
        addition = addition-(bore+access).transform(np.c_[basis, origin])
        wells.append(dict(index=index, axis_local_mm=[sign*7., 0., 0.],
                          inherited_bore_diameter_mm=3.4, added_bore_diameter_mm=3.5,
                          rear_access_diameter_mm=9.,
                          washer_front_z_mm=float(local_washer[:, 2].max()),
                          adapter_bearing_z_mm=bearing_z, nominal_axial_gap_mm=.02))
    adapter = solids['origin:FJ1500/front_support']+addition
    actual, meta = export(adapter, DEST/'integral_saddle_and_grip.stl', origin)
    contacts = []
    for key, body in solids.items():
        if key in ['origin:FJ1500/front_support', 'origin:FJ1500/backing_plate']:
            continue  # Replaced objects, not collision exemptions in the bench assembly.
        volume = float(A.I.intersection_volume(actual, body))
        contacts.append(dict(key=key, volume_mm3=volume))
    loss = float((solids['origin:FJ1500/front_support']-actual).volume())
    wire_volume = float(A.I.intersection_volume(actual, A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), .3)))
    hits = [row for row in contacts if row['volume_mm3'] > .001]
    sources = {key: dict(value) for key, value in trial['sources'].items()
               if key not in ['origin:FJ1500/front_support', 'origin:FJ1500/backing_plate']}
    sources['bench_trial:integral_saddle_and_grip'] = meta
    write_json(DEST/'assembly_check.json', dict(
        scope='Standalone single-muscle zero-pose bench candidate; no bone or distal loading fixture; no drive simulation',
        parent_report=report_file.relative_to(ROOT).as_posix(), parent_sha256=sha(report_file),
        n_instances=len(sources), sources=sources, stage_pass=False, actual_drive_range_deg=None))
    report = dict(scope=__doc__, script_sha256=sha(Path(__file__)),
        current_assembly_sha256=sha(report_file), source_selection_sha256=sha(selection_file),
        trial_source_sha256=sha(trial_file), route_source_sha256=sha(route_file),
        local_basis_columns_source=basis.tolist(), local_origin_source_mm=origin.tolist(),
        original_front_support_loss_mm3=loss, new_adapter=meta,
        export_surface_tolerance_mm=None,
        export_surface_added_signed_mm3=float((actual-adapter).volume()),
        export_surface_removed_signed_mm3=float((adapter-actual).volume()),
        construction_parameters=dict(stem_radius_mm=11., stem_z_interval_mm=[-42., -4.3],
            grip_dimensions_mm=[34., 28., 18.], grip_lower_local_mm=[-17., -14., -56.],
            grip_flat_faces='local X=+-17 or Y=+-14; vice not selected or modelled',
            tpu_relief='subtraction at original and six local-axis shifts of 0.3mm; not a certified uniform clearance',
            active_line_relief_radius_mm=1.2),
        nut_and_tool_wells=wells, contacts=contacts, intersections_over_0_001_mm3=hits,
        active_wire_overlap_mm3=wire_volume, distal_end_constraints_added=False,
        geometric_checks_pass=not hits and loss <= .001 and wire_volume <= .001,
        not_verified=['actual vice and tool fit, screw torque, preload and fastening',
                      'strength, print direction, slicing and material',
                      'active steel terminal gripping and wire friction',
                      'TPU deformation, recovery or arm drive'],
        stage_pass=False, actual_drive_range_deg=None)
    write_json(DEST/'adapter_design.json', report)
    print(json.dumps({key: report[key] for key in ['original_front_support_loss_mm3', 'active_wire_overlap_mm3',
                                                  'intersections_over_0_001_mm3', 'geometric_checks_pass']}), flush=True)
    if not report['geometric_checks_pass']:
        raise ValueError('Bench adapter rejected; preserve evidence and revise geometry')


if __name__ == '__main__':
    main()
