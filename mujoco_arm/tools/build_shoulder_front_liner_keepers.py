"""Positive nominal location/retention for the three preplaced front liners.

Integral TPU nut carrier plus small locating peg; actual rear M2 screw,
spacer and washer. Nut slots are loaded before placing the liner. Original
bone/muscle sources are retained. This does not validate printed fit, preload,
TPU captive-nut pullout strength or complete installation/actuation.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_scapula_clamps as C
import build_shoulder_bench_mount as B
import relieve_shoulder_fixture_for_scapula_insertion as R
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export

SOURCE = R.DEST/'assembly_check.json'
DEST = A.DESIGN/'s1_front_liner_keepers'


def run():
    mf = A.I._manifold()
    DEST.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    sources = dict(parent['sources'])
    solids = {k:A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in sources.items()}
    zones_file = C.DEST/'assembly_check.json'
    zones = json.loads(zones_file.read_text(encoding='utf8'))['zones']
    definitions, modified, new = {}, {}, {}
    for tag, zone in zones.items():
        g = np.array(zone['local_origin_source_mm'])
        F = np.array(zone['basis_columns_source'])
        n = F[:, 2]
        back = zone['outer_front_local_mm']
        key = 'clamp:'+tag+'/front_liner'
        frame = 'fixture:lateral_frame' if tag == 'lateral_upper' else 'fixture:medial_frame'
        old = solids[key]
        local = old.transform(np.c_[F.T, -F.T@g])
        root_rows, additions = [], []
        for name, u, v, radius, end in [('main', 0., 0., 3.3, back+1.),
                                      ('peg', 4.8, 4.8, 1., back-.1)]:
            hits = local.ray_cast([u, v, -20], [u, v, 25])
            z = sorted(float(h.position[2]) for h in hits)
            if len(z) != 2 or z[1]-z[0] < .8:
                raise ValueError('Locator needs an explicit liner attachment interval: '+tag+'/'+name)
            start = z[-1]-.4
            stud = rod([u, v, start], [u, v, end], radius)
            additions.append(stud)
            root_rows.append(dict(name=name, local_uv_mm=[u, v], liner_interval_local_n_mm=z,
                                  start_local_n_mm=start, end_local_n_mm=end, radius_mm=radius,
                                  original_liner_overlap_mm3=A.I.intersection_volume(stud, local)))
            if root_rows[-1]['original_liner_overlap_mm3'] < .1:
                raise ValueError('Locator is not connected to the original lining')
            bore = rod(g+F@np.array([u, v, start-.3]), g+F@np.array([u, v, end+1]), radius+.1)
            solids[frame] -= bore
        # The nut is in the guaranteed rear 3 mm backing, not in the thin
        # skin that contacts the bone. A rearward pillar extension gives
        # 1.65 mm nominal roof above the clearance pocket.
        nut_bottom = back-2.3
        nut_height = 1.6
        nut_radius = 4./np.sqrt(3)
        nut = mf.Manifold.cylinder(nut_height, nut_radius, nut_radius, 6).translate([0, 0, nut_bottom])
        nut -= rod([0, 0, nut_bottom-.2], [0, 0, nut_bottom+nut_height+.2], 1.1)
        pocket = mf.Manifold.cylinder(nut_height+.1, nut_radius+.1, nut_radius+.1, 6).translate([0, 0, nut_bottom-.05])
        entrance = B.box([0, -2.1, nut_bottom-.05], [3.6, 2.1, nut_bottom+nut_height+.05])
        axial_bore = rod([0, 0, back-2.9], [0, 0, back+2.], 1.1)
        studs = additions[0]+additions[1]
        studs -= pocket+entrance+axial_bore
        # Only the additions are bone-trimmed, preserving the existing
        # lining's contact geometry and its independently checked path.
        studs = studs.transform(np.c_[F, g])-solids['bone:scapula']
        revised = old+studs
        solids[key] = revised
        # Open the receiver for the lining's own insertion direction as
        # well. The earlier bone insertion clearance does not imply that
        # a separate liner can be placed through the same contact pocket.
        # A continuous local sweep removes the small lead-in undercuts.
        liner_local = revised.transform(np.c_[F.T, -F.T@g])
        entry = liner_local.minkowski_sum(B.box([-.02, -.02, -20.], [.02, .02, 0.]))
        solids[frame] -= entry.transform(np.c_[F, g])
        # The old exact bearing skin and a newly re-exported receiver can
        # otherwise accumulate a small coplanar material overlap. Make the
        # 0.01 mm normal model fit allowance explicit, preserving the bone
        # side of the liner. It is not a selected printed tolerance.
        solids[frame] -= revised.translate(n*.01)
        modified[key] = dict(added_mm3=float((revised-old).volume()), removed_mm3=float((old-revised).volume()))
        prefix = 'liner_keeper:'+tag+'/'
        # M2x6, a 2.8 mm sleeve, and 0.5 mm washer leave the nominal tip
        # 0.4 mm beyond the nut, still above the bone-side lining surface.
        washer = rod([0, 0, back+2.8], [0, 0, back+3.3], 4.5)-rod([0, 0, back+2.7], [0, 0, back+3.4], 1.1)
        spacer = rod([0, 0, back], [0, 0, back+2.8], 4.5)-rod([0, 0, back-.1], [0, 0, back+2.9], 3.4)
        profile = np.array([[0, back-2.7], [1., back-2.7], [1., back+3.3],
                            [1.9, back+3.3], [1.9, back+5.3], [0, back+5.3]])
        screw = mf.Manifold.revolve(mf.CrossSection([profile]), 64)
        for name, body in dict(nut=nut, screw=screw, washer=washer, spacer=spacer).items():
            new[prefix+name] = body.transform(np.c_[F, g])
        definitions[tag] = dict(liner_key=key, frame_key=frame, source_origin_mm=g.tolist(),
            basis_columns_source=F.tolist(), back_bearing_local_n_mm=back, roots=root_rows,
            main_post_diameter_mm=6.6, main_frame_bore_diameter_mm=6.8,
            secondary_peg_diameter_mm=2., secondary_bore_diameter_mm=2.2,
            nominal_radial_clearance_mm=.1, nut_across_flats_mm=4., nut_height_mm=1.6,
            liner_receiver_normal_model_allowance_mm=.01,
            liner_entry_withdrawal_sweep_mm=20., liner_entry_transverse_model_clearance_mm=.02,
            nut_pocket_radial_clearance_mm=.1, nut_pocket_axial_clearance_each_side_mm=.05,
            nut_bottom_local_n_mm=nut_bottom, nut_entry_direction_local=[1, 0, 0],
            nut_slot_local_bounds_mm=[0, -2.1, nut_bottom-.05, 3.6, 2.1, nut_bottom+1.65],
            nominal_nut_roof_mm=1.65, nominal_screw_length_mm=6., screw_axis_source=(-n).tolist(),
            screw_tip_local_n_mm=back-2.7, nominal_tip_beyond_nut_mm=.4,
            spacer_OD_ID_height_mm=[9., 6.8, 2.8], washer_OD_ID_height_mm=[9., 2.2, .5])
    changed = set(new)|set(modified)|{d['frame_key'] for d in definitions.values()}
    for key, body in {**{k:solids[k] for k in changed if k in solids}, **new}.items():
        positive = [p.volume() for p in body.decompose() if p.volume() > 1e-6]
        if len(positive) != 1:
            (DEST/'rejected_components.json').write_text(json.dumps(dict(key=key, positive_volumes_mm3=positive), indent=2), encoding='utf8')
            raise ValueError('Locator material disconnected: '+key)
        if key in parent['sources']:
            origin = np.array(parent['sources'][key]['translation_source_mm'])
        else:
            bb = np.array(body.bounding_box())
            origin = (bb[:3]+bb[3:])/2
        actual, meta = local_export(body.as_original().simplify(.001), DEST/(key.replace(':', '_').replace('/', '_')+'.stl'), origin)
        if not meta['mesh_checks_pass']:
            (DEST/'rejected_mesh.json').write_text(json.dumps(dict(key=key, meta=meta), indent=2), encoding='utf8')
            raise ValueError('Locator readback fails: '+key)
        solids[key], sources[key] = actual, meta
    changes = {}
    for key in changed & set(parent['sources']):
        row = parent['sources'][key]
        old = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])[0]
        changes[key] = dict(added_mm3=float((solids[key]-old).volume()), removed_mm3=float((old-solids[key]).volume()))
    overlaps = []
    for a, b in itertools.combinations(solids, 2):
        if a not in changed and b not in changed:
            continue
        volume = A.I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            overlaps.append(dict(a=a, b=b, volume_mm3=volume))
    route_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    wires = []
    for row in json.loads(route_file.read_text(encoding='utf8'))['routes']:
        wire = A.D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        volumes = {k:A.I.intersection_volume(wire, s) for k, s in solids.items()}
        wires.append(dict(head=row['head'], volumes_mm3=volumes,
                          hits_over_threshold=[dict(part=k, volume_mm3=v) for k, v in volumes.items() if v > .001]))
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
        parent_report=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
        script_sha256=A.hashed(Path(__file__).resolve()), routes_sha256=A.hashed(route_file),
        generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_scapula_clamps.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py']},
        keeper_definitions=definitions, changed_object_count=len(changed), material_changes=changes,
        changed_object_material_intersections=overlaps, actuation_wire_checks=wires,
        scapula_insertion_relief=dict(direction_source=parent['scapula_insertion_relief']['direction_source'],
                                     maximum_translation_mm=20., inherited_from_parent_report=True),
        intended_body_owners={**parent['intended_body_owners'], **{k:'front_liner_nominal_retainer' for k in new}},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        stage_pass=False, actual_drive_range_deg=None,
        not_verified=['printed fit, captive-nut pullout/rotation, locator bending and actual retaining force',
                      'new liner/nut/spacer/washer/screw insertion and tooling until independently checked',
                      'dense scapula insertion must include all new retainers and be recomputed',
                      'existing physical clamp/material/load and muscle drive limitations still apply'])
    (DEST/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects; changed', len(changed), flush=True)
    print('New overlaps >.001', [r for r in overlaps if r['volume_mm3'] > .001], flush=True)
    print('Active wire failures', sum(bool(w['hits_over_threshold']) for w in wires), flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
