"""Independent two-frame to backboard/bench-interface geometry candidate.

Retains the existing bone-side frame tips. Adds real bores, collars, spacers,
through fasteners and a backboard with feet. The bench patch represents an
assumed externally fixed table region, not a supplied table or solved support.
Bone jaws, tool access, strength and flexible actuation remain unverified.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export

SOURCE = A.DESIGN / 's1_retained_joint_comparison/shaft_collar/station_60/origin_mounts/fit_revision/assembly_check.json'
BASE_SOURCE = A.DESIGN / 's1_port_path_revision/whole_shoulder_static/actual_assembly_check.json'
DEST = SOURCE.parent / 'bench_mount_candidate'
THRESHOLD = .001


def box(lo, hi):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    return A.I._manifold().Manifold.cube((hi-lo).tolist()).translate(lo)


def fasteners(point, axis, bolt_length, far_face, diameter):
    """Plain nominal metric envelopes. Point is the near bearing face.

    Positive axis goes into the joint. Washers contact the bearing faces;
    the thread, compression and preload are deliberately not modeled.
    """
    mf = A.I._manifold()
    p = np.asarray(point, float)
    n = np.asarray(axis, float)
    n /= np.linalg.norm(n)
    u = np.cross(np.eye(3)[np.argmin(abs(n))], n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    basis = np.c_[np.column_stack([u, v, n]), p]
    head_r, head_h, across, nut_h, washer_r = (
        (3.5, 4., 7., 3.2, 4.5) if diameter == 4 else
        (5., 6., 10., 5., 6.))
    hole_r = diameter/2 + (.2 if diameter == 4 else .3)
    wt = .8
    profile = np.array([[0, -wt-head_h], [head_r, -wt-head_h],
                        [head_r, -wt], [diameter/2, -wt],
                        [diameter/2, bolt_length-wt], [0, bolt_length-wt]])
    screw = mf.Manifold.revolve(mf.CrossSection([profile]), 64).transform(basis)
    near = rod(p-n*wt, p, washer_r)-rod(p-n*(wt+.1), p+n*.1, hole_r)
    far = rod(p+n*far_face, p+n*(far_face+wt), washer_r)-rod(
        p+n*(far_face-.1), p+n*(far_face+wt+.1), hole_r)
    nut = mf.Manifold.cylinder(nut_h, across/np.sqrt(3), across/np.sqrt(3), 6)
    nut = nut-mf.Manifold.cylinder(nut_h+.2, diameter/2+.1,
                                 diameter/2+.1, 48).translate([0, 0, -.1])
    nut = nut.translate([0, 0, far_face+wt]).transform(basis)
    projection = bolt_length-wt-(far_face+wt+nut_h)
    if projection <= 0:
        raise ValueError('Bolt does not emerge through full nut envelope')
    return dict(screw=screw, near_washer=near, far_washer=far, nut=nut), dict(
        nominal_diameter_mm=diameter, nominal_length_mm=bolt_length,
        axis_source=n.tolist(), near_bearing_source_mm=p.tolist(),
        far_bearing_source_mm=(p+n*far_face).tolist(),
        nominal_tip_beyond_nut_mm=projection,
        thread_and_preload_verified=False)


def run():
    mf = A.I._manifold()
    DEST.mkdir(parents=True, exist_ok=True)
    parent = json.loads(SOURCE.read_text(encoding='utf8'))
    old = json.loads(BASE_SOURCE.read_text(encoding='utf8'))
    solids = {k: A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
              for k, r in parent['sources'].items()}
    sources = dict(parent['sources'])
    # The standing plate is behind the shoulder (+Y); feet project further
    # behind it so the hanging arm is not placed through the tabletop patch.
    board = box([-180, 30, 1222], [45, 36, 1370])
    bench = box([-185, 44, 1198], [50, 134, 1216])
    for x in [-155., 5.]:
        board = board + box([x-15, 35, 1216], [x+15, 110, 1224])
        # Put two braces either side of the M6 screw/washer, not through
        # their installation axis. Keep both inside the 30 mm foot width.
        for side in [-1, 1]:
            bx = x+side*10.5
            board = board + rod([bx, 34, 1260], [bx, 94, 1220], 4.)
    new, mounts, contacts, bores = {}, {}, [], []
    original_frames = {}
    for name in ['medial', 'lateral']:
        key = 'fixture:'+name+'_frame'
        original_frames[key] = solids[key]
        r = old['sources']['fixture:'+name+'_base']
        base = A.load_checked(A.ROOT/r['file'], r['translation_source_mm'], r['sha256'])[0]
        if name == 'lateral':
            # This displacement belongs to the existing retained frame;
            # the frame/bone points themselves are not relocated here.
            base = base.translate([-9., 0, 0])
        bb = np.array(base.bounding_box())
        z = (bb[2]+bb[5])/2
        # Offset lateral holes away from the oblique post's existing root.
        offsets = [12., 30.] if name == 'medial' else [22., 36.]
        frame = solids[key]
        records = []
        for i, dx in enumerate(offsets):
            p = np.array([bb[0]+dx, bb[1], z])
            q = np.array([p[0], bb[4], z])
            # Integrated R6 bearing bosses increase the local Z section.
            frame = frame + rod(p, q, 6.)
            drill = rod(p-[0, 2, 0], [p[0], 45, z], 2.2)
            frame = frame-drill
            board = board-drill
            spacer = rod(q, [q[0], 30, z], 6.)-drill
            label = f'bench:{name}_{i}'
            new[label+'/spacer'] = spacer
            hardware, spec = fasteners(p, [0, 1, 0], 60., 36.-p[1], 4)
            new.update({label+'/'+k: s for k, s in hardware.items()})
            records.append(dict(index=i, centre_source_mm=p.tolist(),
                                spacer_length_mm=float(30.-q[1]), fastener=spec))
            contacts.extend([(key, label+'/spacer', [0, 1, 0]),
                             (label+'/spacer', 'bench:backboard', [0, 1, 0]),
                             (label+'/near_washer', key, [0, 1, 0]),
                             (label+'/far_washer', 'bench:backboard', [0, -1, 0])])
            bores.append(dict(name=label, start_source_mm=(p-[0, 2, 0]).tolist(),
                              end_source_mm=[p[0], 45, z], shaft_radius_mm=2.,
                              parts=[key, label+'/spacer', 'bench:backboard']))
        solids[key] = frame
        mounts[name] = dict(base_bounds_source_mm=bb.tolist(), bores=records)
    for i, x in enumerate([-155., 5.]):
        p = np.array([x, 90., 1224.])
        bore = rod([x, 90, 1184], [x, 90, 1235], 3.3)
        board = board-bore
        bench = bench-bore
        label = f'bench:foot_{i}'
        hardware, spec = fasteners(p, [0, 0, -1], 35., 26., 6)
        new.update({label+'/'+k: s for k, s in hardware.items()})
        mounts[f'foot_{i}'] = spec
        contacts.extend([('bench:backboard', 'bench:table_interface_patch', [0, 0, -1]),
                         (label+'/near_washer', 'bench:backboard', [0, 0, -1]),
                         (label+'/far_washer', 'bench:table_interface_patch', [0, 0, 1])])
        bores.append(dict(name=label, start_source_mm=[x, 90, 1184],
                          end_source_mm=[x, 90, 1235], shaft_radius_mm=3.,
                          parts=['bench:backboard', 'bench:table_interface_patch']))
    new.update({'bench:backboard': board, 'bench:table_interface_patch': bench})
    changes = {}
    for key, solid in {**{k: solids[k] for k in original_frames}, **new}.items():
        bb = np.array(solid.bounding_box())
        origin = (bb[:3]+bb[3:])/2
        path = DEST/(key.replace(':', '_').replace('/', '_')+'.stl')
        restored, meta = local_export(solid.as_original().simplify(.001), path, origin)
        if not meta['mesh_checks_pass']:
            (DEST/'rejected_mesh.json').write_text(json.dumps(dict(key=key, meta=meta), indent=2), encoding='utf8')
            raise ValueError('Bench candidate mesh fails: '+key)
        solids[key], sources[key] = restored, meta
        if key in original_frames:
            changes[key] = dict(added_mm3=float((restored-original_frames[key]).volume()),
                                removed_mm3=float((original_frames[key]-restored).volume()))
    # Check every pair containing any changed/new object; unchanged parent
    # pairs are carried in a separate field with its exact source hash.
    changed = set(new) | set(original_frames)
    intersections = []
    for a, b in itertools.combinations(solids, 2):
        if a not in changed and b not in changed:
            continue
        vol = A.I.intersection_volume(solids[a], solids[b])
        if vol > 1e-9:
            intersections.append(dict(a=a, b=b, volume_mm3=vol,
                                      flexible=any(k.startswith('tpu:') or k.startswith('collar:liner') for k in [a, b])))
    route_path = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes = json.loads(route_path.read_text(encoding='utf8'))['routes']
    wire_checks = []
    for route in routes:
        wire = A.D.outer_envelope(np.array(route['full_wire_path_source_mm']), .3)
        hits = [dict(part=k, volume_mm3=A.I.intersection_volume(wire, s)) for k, s in solids.items()]
        wire_checks.append(dict(head=route['head'], intersections=hits,
                                hits_over_threshold=[h for h in hits if h['volume_mm3'] > THRESHOLD]))
    contact_checks = []
    for a, b, axis in contacts:
        direction = np.array(axis, float)
        contact_checks.append(dict(a=a, b=b, gap_capped_1mm=A.I.gap(solids[a], solids[b], 1.),
                                   nominal_intersection_mm3=A.I.intersection_volume(solids[a], solids[b]),
                                   probe_translation_source_mm=(direction*.05).tolist(),
                                   intersection_after_0_05mm_mm3=A.I.intersection_volume(solids[a].translate(direction*.05), solids[b])))
    bore_checks = []
    for bore in bores:
        probe = rod(bore['start_source_mm'], bore['end_source_mm'], bore['shaft_radius_mm'])
        bore_checks.append({**bore, 'shaft_material_intersections_mm3': {
            k:A.I.intersection_volume(probe, solids[k]) for k in bore['parts']}})
    owners = dict(parent['intended_body_owners'])
    owners.update({k: 'nominal_bench_fastened' for k in new})
    owners['bench:table_interface_patch'] = 'assumed_external_fixed_bench_region'
    report = dict(scope=__doc__, n_instances=len(solids), sources=sources,
                  parent_report=SOURCE.relative_to(A.ROOT).as_posix(), parent_sha256=A.hashed(SOURCE),
                  original_base_source_sha256=A.hashed(BASE_SOURCE), script_sha256=A.hashed(Path(__file__).resolve()),
                  route_sha256=A.hashed(route_path), base_mounts=mounts, frame_modifications=changes,
                  intended_body_owners=owners,
                  free_terminal_constraints_required=parent['free_terminal_constraints_required'],
                  changed_object_count=len(changed), changed_object_material_intersections=intersections,
                  unchanged_parent_material_intersections=parent['material_intersections'],
                  actuation_wire_checks=wire_checks, finite_contact_checks=contact_checks,
                  shaft_bore_checks=bore_checks,
                  threshold_mm3=THRESHOLD,
                  static_new_pairs_below_threshold=all(r['volume_mm3'] <= THRESHOLD for r in intersections),
                  stage_pass=False, actual_drive_range_deg=None,
                  not_verified=['bone contact shoes, scapula/clavicle jaws and their fastening still incomplete',
                                'bench patch is assumed support, not supplied table dimensions or anchorage',
                                'bolt threads, preload, friction, strength, printability and tool access',
                                'assembly insertion/withdrawal paths and rigidity under load',
                                'same unresolved compliant contact, channel-wall differences and free-terminal dynamics as parent',
                                'no new rigid or flexible motion range or physical drive result'])
    (DEST/'assembly_check.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Readback', len(solids), 'objects; changed', len(changed), flush=True)
    print('New overlaps > threshold:', [r for r in intersections if r['volume_mm3'] > THRESHOLD], flush=True)
    print('Active wire failures:', sum(bool(w['hits_over_threshold']) for w in wire_checks), flush=True)
    print('Contact probes absent:', [r for r in contact_checks if r['intersection_after_0_05mm_mm3'] <= THRESHOLD], flush=True)
    print('Blocked bores:', [r for r in bore_checks if max(r['shaft_material_intersections_mm3'].values()) > THRESHOLD], flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
