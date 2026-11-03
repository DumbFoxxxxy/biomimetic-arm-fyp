"""Split clavicle clamp and removable support to the existing S1 backboard.

Uses section-centroid regression to locate a local fixture frame. Keeps all
bone/muscle meshes unchanged. Nominal shells, liners and through fasteners
are actual meshes, but preload, friction, strength and printability are not
validated. The support is a test-bench candidate, not a mobile SC joint.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_scapula_clamps as C
import build_shoulder_bench_mount as B
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export

SOURCE = C.DEST/'assembly_check.json'
DEST = C.DEST/'clavicle_mount_candidate'


def centroid(section):
    area, moment = 0., np.zeros(2)
    for p in section.to_polygons():
        q = np.roll(p, -1, axis=0)
        cross = p[:, 0]*q[:, 1]-q[:, 0]*p[:, 1]
        area += cross.sum()/2
        moment += ((p+q)*cross[:, None]).sum(axis=0)/6
    if area <= .01:
        raise ValueError('Empty or inconsistently oriented clavicle section')
    return moment/area


def local_frame(bone, station):
    joints = A.frames.joint_centres('source')
    z = A.D.unit(joints['AC']-joints['SC'])
    y = A.D.unit(np.array([0., 1, 0])-z*z[1])
    x = np.cross(y, z)
    F = np.column_stack([x, y, z])
    local = bone.transform(np.c_[F.T, -F.T@joints['SC']])
    sample_z = np.arange(station-6, station+6.01, 1.)
    centres = np.array([centroid(local.slice(float(q))) for q in sample_z])
    fit = np.array([np.polyfit(sample_z-station, centres[:, k], 1) for k in [0, 1]])
    prediction = np.c_[sample_z-station, np.ones(len(sample_z))]@fit.T
    origin = joints['SC']+F@np.r_[fit[:, 1], station]
    tangent = A.D.unit(F@np.r_[fit[:, 0], 1.])
    posterior = A.D.unit(np.array([0., 1, 0])-tangent*tangent[1])
    cross = np.cross(posterior, tangent)
    frame = np.column_stack([cross, posterior, tangent])
    return origin, frame, dict(
        station_on_SC_AC_reference_axis_mm=station,
        method='area-weighted section centroids; linear regression over station +/-6 mm, 1 mm pitch',
        SC_source_mm=joints['SC'].tolist(), AC_source_mm=joints['AC'].tolist(),
        reference_basis_columns_source=F.tolist(), sample_stations_mm=sample_z.tolist(),
        section_centroids_reference_xy_mm=centres.tolist(), fit_coefficients_xy=fit.tolist(),
        centroid_fit_xy_RMS_mm=float(np.sqrt(np.mean((centres-prediction)**2))),
        source_origin_mm=origin.tolist(), basis_columns_source=frame.tolist(),
        not_an_anatomical_joint_axis=True)


def run(station, dest):
    mf = A.I._manifold()
    dest.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    sources = dict(parent['sources'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in sources.items()}
    bone = solids['bone:clavicle']
    g, F, frame_meta = local_frame(bone, station)
    local_bone = bone.transform(np.c_[F.T, -F.T@g])
    profile = mf.CrossSection.batch_boolean(
        [local_bone.slice(float(h)) for h in np.arange(-6., 6.01, .5)], mf.OpType.Add).hull()
    inner = profile.offset(1., circular_segments=48)
    rigid_inner = inner.offset(.02, circular_segments=48)
    outer = rigid_inner.offset(3., circular_segments=48)
    bb = np.array(outer.bounds())
    cx, cy = (bb[:2]+bb[2:])/2
    cavity = rigid_inner.extrude(16.).translate([0, 0, -8.])
    shell = outer.extrude(12.).translate([0, 0, -6.])-cavity
    liner_blank = inner.extrude(11.6).translate([0, 0, -5.8])
    local_stock = local_bone ^ B.box([-40, -40, -6.1], [40, 40, 6.1])
    parts, definitions, jaws = {}, {}, {}
    bolt_x = [bb[0]-3.5, bb[2]+3.5]
    for sign, label in [(1, 'posterior'), (-1, 'anterior')]:
        # Open each half toward the opposite half. A merely conformal
        # closed bone pocket can contain undercuts, even inside an outer
        # convex clamp profile. Keep this a geometric opening operation,
        # never a fictitious rigidly attached/deforming soft liner.
        swept = B.box([-.02, -20 if sign > 0 else 0, -.02],
                      [.02, 0 if sign > 0 else 20, .02])
        bone_opening = local_stock.minkowski_sum(swept)
        rigid_opening = cavity.minkowski_sum(swept)
        jaw = shell
        for bx in bolt_x:
            jaw += rod([bx, cy+sign*.5, 0], [bx, cy+sign*6, 0], 3.7)
            jaw += rod([cx, cy+sign*3., 0], [bx, cy+sign*3., 0], 2.5)
        clip = B.box([-100, cy+.5 if sign > 0 else -100, -10],
                     [100, 100 if sign > 0 else cy-.5, 10])
        jaw = (jaw-rigid_opening)^clip
        for bx in bolt_x:
            jaw -= rod([bx, cy-20, 0], [bx, cy+20, 0], 1.7)
            # The rounded ear bridges must not project into the real
            # washer seating planes at +/-6 mm.
            jaw -= rod([bx, cy+sign*6, 0], [bx, cy+sign*20, 0], 3.75)
        jaws[label] = jaw
        parts['liner_'+label] = ((liner_blank-bone_opening)^clip).transform(np.c_[F, g])
    for i, bx in enumerate(bolt_x):
        near = g+F@np.array([bx, cy+6., 0])
        hardware = C.fasteners(near, -F[:, 1], 20., 12.)
        parts.update({f'M3_{k}_{i}':s for k, s in hardware.items()})
        definitions[f'M3_{i}'] = dict(near_bearing_source_mm=near.tolist(), direction_source=(-F[:, 1]).tolist(),
                                      nominal_length_mm=20., required_withdrawal_mm=20.1,
                                      bore_axis_source_mm=(g+F@np.array([bx, cy, 0])).tolist())
    # Posterior side is supported toward the backboard (+source Y).
    # The central member and two lower diagonal braces form one support
    # part; materials/load capacity and manufacture remain to be chosen.
    start = g+F@np.array([cx, bb[3]-1.5, 0])
    q = np.array([start[0], 24., start[2]])
    if q[1]-start[1] < 20:
        raise ValueError('Unexpected clavicle/backboard orientation')
    support = rod(start-[0, 1, 0], q+[0, 1, 0], 5.)
    branch = start+np.array([0, 15., 0])
    for sign in [-1, 1]:
        support += rod(branch, q+np.array([sign*17, 1, -10]), 3.5)
    flange = B.box([q[0]-23, 24, q[2]-14], [q[0]+23, 30, q[2]+14])
    fixed = jaws['posterior'].transform(np.c_[F, g])+support+flange
    posterior_opening = cavity.minkowski_sum(B.box([-.02, -20, -.02], [.02, 0, .02]))
    fixed -= posterior_opening.transform(np.c_[F, g])
    parts['moving_jaw'] = jaws['anterior'].transform(np.c_[F, g])
    board = solids['bench:backboard']
    for i, dx in enumerate([-13., 13.]):
        near = q+np.array([dx, 0, 0])
        bore = rod(near-[0, 2, 0], near+[0, 22, 0], 2.2)
        fixed -= bore
        board -= bore
        hardware, dimensions = B.fasteners(near, [0, 1, 0], 20., 12., 4)
        parts.update({f'M4_{k}_{i}':s for k, s in hardware.items()})
        definitions[f'M4_{i}'] = dimensions
    # Real bores remain open through the supporting material as well.
    for bx in bolt_x:
        p = g+F@np.array([bx, cy, 0])
        fixed -= rod(p-F[:, 1]*20, p+F[:, 1]*20, 1.7)
    parts['fixed_support'] = fixed
    generated = {'clavicle:'+k:s for k, s in parts.items()}
    generated['bench:backboard'] = board
    for key, solid in generated.items():
        positive = [p.volume() for p in solid.decompose() if p.volume() > 1e-6]
        if len(positive) != 1:
            (dest/'rejected_components.json').write_text(json.dumps(dict(key=key, positive_volumes_mm3=positive), indent=2), encoding='utf8')
            raise ValueError('Disconnected clavicle fixture: '+key)
        b = np.array(solid.bounding_box())
        origin = (b[:3]+b[3:])/2
        actual, meta = local_export(solid.as_original().simplify(.001), dest/(key.replace(':', '_')+'.stl'), origin)
        if not meta['mesh_checks_pass']:
            (dest/'rejected_mesh.json').write_text(json.dumps(dict(key=key, meta=meta), indent=2), encoding='utf8')
            raise ValueError('Clavicle fixture readback fails: '+key)
        solids[key], sources[key] = actual, meta
    contacts = []
    for a, b in itertools.combinations(solids, 2):
        if a not in generated and b not in generated:
            continue
        vol = A.I.intersection_volume(solids[a], solids[b])
        if vol > 1e-9:
            contacts.append(dict(a=a, b=b, volume_mm3=vol))
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    wire_checks = []
    for row in json.loads(route_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        hits = [dict(part=k, volume_mm3=A.I.intersection_volume(wire, s)) for k, s in solids.items()]
        wire_checks.append(dict(head=row['head'], intersections=hits,
                                hits_over_threshold=[h for h in hits if h['volume_mm3'] > .001]))
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
        parent_report=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
        script_sha256=A.hashed(Path(__file__).resolve()), routes_sha256=A.hashed(route_file),
        generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [A.ROOT/'mujoco_arm/tools/build_shoulder_scapula_clamps.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py']},
        clavicle_fixture_frame=frame_meta,
        dimensions=dict(axial_length_mm=12., liner_axial_length_mm=11.6, profile_sample_step_mm=.5,
                        liner_envelope_offset_mm=1., rigid_shell_wall_mm=3., split_gap_mm=1.,
                        liner_outer_to_jaw_profile_clearance_mm=.02,
                        half_opening_sweep_mm=20., half_opening_transverse_clearance_mm=.02,
                        M3_ear_bearing_half_spacing_mm=6., M3_ear_bore_diameter_mm=3.4,
                        main_support_radius_mm=5., brace_radius_mm=3.5,
                        flange_source_min_mm=[q[0]-23,24,q[2]-14], flange_source_max_mm=[q[0]+23,30,q[2]+14]),
        clamp_cross_section_bounds_local_mm=bb.tolist(), local_bolt_x_mm=list(bolt_x),
        fasteners=definitions, main_support_start_source_mm=start.tolist(), flange_centre_near_source_mm=q.tolist(),
        main_member_length_mm=float(q[1]-start[1]+2),
        changed_object_count=len(generated), changed_object_material_intersections=contacts,
        actuation_wire_checks=wire_checks,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**parent['intended_body_owners'],
                              **{k:'nominal_clavicle_bench_clamp' for k in generated if k != 'bench:backboard'}},
        static_new_pairs_below_threshold=all(r['volume_mm3'] <= .001 for r in contacts),
        stage_pass=False, actual_drive_range_deg=None,
        not_verified=['clamp preload, friction, liner compression, support bending and bone/fastener strength',
                      'material, print/metal fabrication, tolerances and actual liner holding',
                      'bone/liner/jaw insertion paths, nut/washer and tool access until independently checked',
                      'same unresolved muscle soft contact, channel-wall residuals and coupled free-terminal actuation'])
    (dest/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects; changed', len(generated), flush=True)
    print('New overlaps >.001', [r for r in contacts if r['volume_mm3'] > .001], flush=True)
    print('Active wire failures', sum(bool(r['hits_over_threshold']) for r in wire_checks), flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser = argparse.ArgumentParser()
    parser.add_argument('--station', type=float, default=15.)
    parser.add_argument('--out', type=Path, default=DEST)
    args = parser.parse_args()
    run(args.station, args.out)
