"""Reusable METAL screw-stop candidate for the current FJ1500 active wire.

Retains the existing counterbore and cap. Radial M2x0.4 screw is a smooth
clearance envelope with a flat point, not an explicitly threaded solid.
First-touch cable state does not establish clamping force or grip.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'mujoco_arm/tools'), str(ROOT)]
import audit_shoulder_retained_joint_fit as A
import design_shoulder_tendon_link as L
import mesh_geom
from build_s1_infraspinatus_bench import export

DEST = ROOT/'prototype/s1_infraspinatus_trial/screw_stop_v1'


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    assembly_file = ROOT/selection['assembly_report']
    if sha(assembly_file) != selection['assembly_sha256']:
        raise ValueError('Stale selected assembly')
    assembly = read(assembly_file)
    route_file = ROOT/read(DEST.parent/'manifest.json')['active_route_source']
    route = next(row for row in read(route_file)['routes'] if row['head'] == 'FJ1500')
    spec = route['distal_terminal_override']
    origin = np.array(spec['inner_face_source_mm'])
    normal = A.D.unit(spec['normal_source'])
    tangent = A.D.unit(spec['tangent_source'])
    basis = np.stack([-tangent, -np.cross(normal, tangent), normal], axis=1)
    assert np.allclose(basis.T@basis, np.eye(3)) and np.linalg.det(basis) > .999999
    mf = A.I._manifold()
    def cylinder(z0, z1, radius):
        return mf.Manifold.cylinder(z1-z0, radius, radius, 96).translate([0., 0., z0])
    def radial_z_to_x(solid, offset=(0., 0., 5.05)):
        matrix = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
        return solid.transform(np.c_[matrix, np.array(offset)])
    # The through-hole and radial thread have manufacturing callouts in the
    # drawing. A 2.1mm cylinder represents clearance, NOT a tap-drill size.
    blank = cylinder(.45, 3.05, 1.25)+cylinder(3.05, 7.05, 1.8)
    local_body = blank-cylinder(.2, 7.3, .4)
    radial_port = radial_z_to_x(mf.Manifold.cylinder(3., 1.05, 1.05, 64))
    local_body = local_body-radial_port
    # Catalog-type first-touch envelope: M2x3, 0.9mm flat point, 0.9 AF
    # hex socket, 0.8mm socket depth. Thread crests/flanks are not simulated.
    local_screw = (mf.Manifold.cylinder(.3, .45, 1., 64).translate([0., 0., .3])+
                   mf.Manifold.cylinder(2.7, 1., 1., 64).translate([0., 0., .6]))
    socket = mf.Manifold.cylinder(1., .9/np.sqrt(3), .9/np.sqrt(3), 6).translate([0., 0., 2.5])
    local_screw = radial_z_to_x(local_screw-socket)
    affine = np.c_[basis, origin]
    body, body_meta = export(local_body.transform(affine), DEST/'METAL_stop_body_envelope.stl', origin)
    screw, screw_meta = export(local_screw.transform(affine), DEST/'M2x3_flat_point_NOT_PRINT.stl', origin)
    # Extend only the tail beyond the existing distal outlet. No internal
    # route, anchor, cap, TPU or proximal wire geometry changes.
    path = np.array(route['full_wire_path_source_mm'])
    extension_end = np.array(spec['inner_face_source_mm'])+normal*8.05
    extended_path = np.vstack([path, extension_end])
    wire, wire_meta = export(L.wire_tube(extended_path), DEST/'active_wire_extended_NOT_PRINT.stl', origin)
    sources = {}
    solids = {}
    old_stop_key = 'stop:FJ1500/distal_stop_blank'
    for key, source in assembly['sources'].items():
        if key == old_stop_key:
            continue  # A replacement, not an ignored colliding pair.
        solids[key] = A.load_checked(ROOT/source['file'], source['translation_source_mm'], source['sha256'])[0]
        sources[key] = source
    conflicts = []
    for new_key, new_solid in [('stop_body', body), ('set_screw', screw), ('extended_active_wire', wire)]:
        for old_key, old in solids.items():
            volume = float(A.I.intersection_volume(new_solid, old))
            if volume > .001:
                conflicts.append(dict(a=new_key, b=old_key, volume_mm3=volume))
    internal = [dict(a='stop_body', b='set_screw', volume_mm3=float(A.I.intersection_volume(body, screw))),
                dict(a='stop_body', b='active_wire', volume_mm3=float(A.I.intersection_volume(body, wire))),
                dict(a='set_screw', b='active_wire', volume_mm3=float(A.I.intersection_volume(screw, wire)))]
    conflicts.extend(row for row in internal if row['volume_mm3'] > .001)
    sources.update({'stop_screw:FJ1500/body': body_meta,
                    'stop_screw:FJ1500/M2x3': screw_meta,
                    'wire_active:FJ1500/extended_reference': wire_meta})
    full = dict(scope='Full current shoulder plus screw-stop first-touch candidate and active wire; static only',
                sources=sources, n_instances=len(sources), stage_pass=False, actual_drive_range_deg=None)
    write(DEST/'full_shoulder_check.json', full)
    bench_file = DEST.parent/'bench_adapter_v1/assembly_check.json'
    bench = read(bench_file)
    bench_sources = {key: source for key, source in bench['sources'].items() if key != old_stop_key}
    bench_sources.update({key: sources[key] for key in sources if key.startswith(('stop_screw:', 'wire_active:'))})
    write(DEST/'bench_check.json', dict(scope='Bench plus screw-stop and extended active wire, first-touch only',
        sources=bench_sources, n_instances=len(bench_sources), stage_pass=False, actual_drive_range_deg=None))
    design = dict(scope=__doc__, source_assembly=selection['assembly_report'], source_assembly_sha256=sha(assembly_file),
        script_sha256=sha(Path(__file__)), source_route=route_file.relative_to(ROOT).as_posix(), source_route_sha256=sha(route_file),
        bench_source_sha256=sha(bench_file), local_origin_source_mm=origin.tolist(), local_basis_columns_source=basis.tolist(),
        nominal_dimensions=dict(stub_diameter_mm=2.5, stub_z_interval_mm=[.45, 3.05],
            head_diameter_mm=3.6, head_z_interval_mm=[3.05, 7.05], through_bore_mm=.8,
            radial_thread='M2x0.4; envelope only, tool and thread engagement unverified',
            radial_port_envelope_mm=2.1, radial_thread_axis_local_z_mm=5.05,
            radial_thread_axis_local_direction=[1., 0., 0.], set_screw_nominal_length_mm=3.,
            nominal_screw_tip_x_mm=.3, first_touch=True, minimum_nominal_screw_body_engagement_mm=1.2,
            flat_point_nominal_diameter_mm=.9, socket_nominal_AF_mm=.9, socket_nominal_depth_mm=.8,
            axial_seating_gap_mm=.05),
        alternatives=dict(rejected_head_diameter_mm=4., observed_cap_overlap_mm3=.24126474283949728,
                          selected_head_diameter_mm=3.6, comparison='neutral source geometry only; no strength ranking'),
        components={'body': body_meta, 'set_screw': screw_meta, 'active_wire': wire_meta},
        extended_active_wire_path_source_mm=extended_path.tolist(),
        added_tail_length_mm=float(np.linalg.norm(extension_end-path[-1])),
        internal_route_and_anchors_unchanged=True, static_conflicts_over_0_001_mm3=conflicts,
        static_internal_intersections=internal,
        screw_catalog_sources=[
            'https://www.pts-uk.com/products/socket-screws/socket-set-screw-flat-point/metric-a2/a9130203',
            'https://www.del-screw.co.jp/en/imo/hirasaki/fe.html'],
        static_geometry_pass=not conflicts, grip_force_N=None, tightening_torque_Nm=None,
        not_verified=['custom metal body manufacture/material/thread fit and physical tool access',
                      'clamping torque, grip, steel strand damage, loosening and repeatability',
                      'printed cap bearing, loaded TPU motion or full-arm actuation'],
        stage_pass=False, actual_drive_range_deg=None)
    write(DEST/'stop_design.json', design)
    print(json.dumps(dict(static_geometry_pass=design['static_geometry_pass'], conflicts=conflicts,
                          full_scene_parts=len(sources), bench_scene_parts=len(bench_sources))), flush=True)
    if conflicts:
        raise ValueError('Stop candidate conflicts with actual source geometry; see saved report')


if __name__ == '__main__':
    main()
