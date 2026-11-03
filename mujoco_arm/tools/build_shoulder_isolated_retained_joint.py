"""Independent GH-only comparison candidate, preserving prior files.

Apply the G1 GH geometry delta to the anatomical scapula, not the older
AC-modified scapula. Rebuild two defective nominal screw envelopes with a
single revolved profile. Export/read back every new part before checking the
complete available zero-pose assembly and six existing wire paths.
"""
import argparse
import itertools
import json
import shutil
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import frames
import interference as I
import mesh_geom
import design_shoulder_port_paths as D
from build_shoulder_tendon_segments import rod

DEST = A.DESIGN / 's1_retained_joint_comparison/GH_only_candidate'


def screw_envelope(site, axis, depth, clearance=0.):
    """Nominal M3x12 envelope; no thread or fastener strength claim."""
    z = np.asarray(axis, float)
    z = z / np.linalg.norm(z)
    helper = np.eye(3)[np.argmin(abs(z))]
    x = np.cross(helper, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    head = depth + .75
    # One closed cross-section avoids coincident cylinder caps at the step.
    profile = [[0, head-12-clearance], [1.5+clearance, head-12-clearance],
               [1.5+clearance, head-clearance], [2.75+clearance, head-clearance],
               [2.75+clearance, head+3+clearance], [0, head+3+clearance]]
    mf = I._manifold()
    shape = mf.Manifold.revolve(mf.CrossSection([np.array(profile)]), 64)
    return shape.transform(np.c_[np.column_stack([x, y, z]), np.asarray(site)])


def export_checked(solid, path):
    display = solid.translate(-frames.SHIFT)
    v, f = mesh_geom.manifold_mesh_to_arrays(display.to_mesh64())
    mesh_geom.write_stl(str(path), v, f)
    result, tri, meta = A.load_checked(path, frames.SHIFT)
    meta['positive_components'] = len([c for c in result.decompose() if c.volume() > 1e-6])
    if meta['positive_components'] != 1:
        raise ValueError('Disconnected exported part: '+str(path))
    meta['export_volume_change_mm3'] = float(result.volume()-solid.volume())
    if abs(meta['export_volume_change_mm3']) > max(.01, solid.volume()*1e-4):
        raise ValueError('Excessive STL volume change: '+str(path))
    return result, meta


def run(dest):
    dest.mkdir(parents=True, exist_ok=True)
    joint_file = A.G1/'assembly.json'
    assembly_file = A.CURRENT/'actual_assembly_check.json'
    routes_file = A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    joint = json.loads(joint_file.read_text(encoding='utf8'))
    assembly = json.loads(assembly_file.read_text(encoding='utf8'))
    routes = json.loads(routes_file.read_text(encoding='utf8'))['routes']
    solids, sources = {}, {}
    for key, row in assembly['sources'].items():
        s, _, meta = A.load_checked(A.ROOT/row['file'], row['translation_source_mm'], row['sha256'])
        solids[key] = s
        sources[key] = meta
    oldscap = solids.pop('bone:scapula')
    oldhum = solids.pop('bone:humerus')
    sources.pop('bone:scapula')
    sources.pop('bone:humerus')
    ac_path = A.G1.parent/'scapula_AC_socket.obj'
    base, _, base_meta = A.load_checked(ac_path, frames.SHIFT)
    engineered_path = A.G1/'scapula_AC_GH_sockets.obj'
    engineered, _, eng_meta = A.load_checked(engineered_path, frames.SHIFT)
    added = engineered-base
    removed = base-engineered
    raw = oldscap-removed+added
    # The delta contains nearly coincident copies of the original bone skin.
    # Simplify their sub-tolerance edges before float32 STL conversion; measure
    # all material changes and still require a closed, single-component readback.
    prepared = raw.as_original().simplify(.01)
    vertices, faces = mesh_geom.manifold_mesh_to_arrays(prepared.to_mesh64())
    vertices, faces = mesh_geom.weld(vertices, faces)
    edges = np.sort(np.concatenate([faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    fourway = unique[counts > 2]
    seam_repairs = []
    if len(fourway):
        if len(fourway) > 2 or np.any(counts[counts > 2] != 4):
            raise ValueError('Unexpected seam topology; explicit redesign required')
        failure = dest/'iterations/unrounded_seam'
        failure.mkdir(parents=True, exist_ok=True)
        fv, ff = mesh_geom.manifold_mesh_to_arrays(prepared.translate(-frames.SHIFT).to_mesh64())
        mesh_geom.write_stl(str(failure/'scapula_GH_only.stl'),fv,ff)
        for edge in fourway:
            p, q = vertices[edge]
            direction = (q-p)/np.linalg.norm(q-p)
            if np.linalg.norm(q-p) > 4:
                raise ValueError('Seam exceeds local rounding scope')
            # Add an explicit R0.1 seam fillet instead of deleting material or
            # silently treating four faces on one edge as a printable solid.
            rounding = rod(p-.1*direction, q+.1*direction, .1)
            before = float(prepared.volume())
            prepared = prepared + rounding
            seam_repairs.append(dict(edge_source_mm=[p.tolist(),q.tolist()],
                radius_mm=.1, extension_mm=.1,
                added_mm3=float(prepared.volume())-before))
        prepared = prepared.as_original().simplify(.001)
    isolated, meta = export_checked(prepared, dest/'scapula_GH_only.stl')
    meta['simplification_tolerance_mm'] = .01
    meta['simplification_added_mm3'] = float((prepared-raw).volume())
    meta['simplification_removed_mm3'] = float((raw-prepared).volume())
    meta['local_added_seam_rounding'] = seam_repairs
    solids['bone:scapula_GH_only'] = isolated
    sources['bone:scapula_GH_only'] = meta
    for key in joint['components']:
        if key == 'scapula_AC_GH_sockets':
            continue
        if key.startswith('GH_bolt_'):
            index = int(key.rsplit('_', 1)[1])
            site = np.array(joint['bolt_sites_display_mm'][index])+frames.SHIFT
            s = screw_envelope(site, joint['neck_axis'], joint['cap_depth_mm'])
            s, meta = export_checked(s, dest/(key+'_profile.stl'))
        else:
            s, _, meta = A.load_checked(A.G1/(key+'.obj'), frames.SHIFT)
        solids[key] = s
        sources[key] = meta
    bolt_reliefs = []
    for index in range(2):
        key = 'GH_bolt_'+str(index)
        overlap = I.intersection_volume(isolated, solids[key])
        if overlap <= 1e-6:
            continue
        archive = dest/'iterations/before_screw_head_relief'
        archive.mkdir(parents=True, exist_ok=True)
        if not (archive/'scapula_GH_only.stl').exists():
            shutil.copyfile(dest/'scapula_GH_only.stl',archive/'scapula_GH_only.stl')
        site = np.array(joint['bolt_sites_display_mm'][index])+frames.SHIFT
        clearance = screw_envelope(site,joint['neck_axis'],joint['cap_depth_mm'],.15)
        revised = isolated-clearance
        before = float(isolated.volume())
        isolated, revised_meta = export_checked(revised,dest/'scapula_GH_only.stl')
        bolt_reliefs.append(dict(screw=key,overlap_before_mm3=overlap,
            nominal_clearance_mm=.15,removed_mm3=before-float(isolated.volume()),
            overlap_after_mm3=I.intersection_volume(isolated,solids[key])))
        sources['bone:scapula_GH_only'].update(revised_meta)
    solids['bone:scapula_GH_only'] = isolated
    # Keep the canonical bone keys so the existing actual-mesh motion checker
    # assigns the engineered head to the humerus rigid body, not the fixture.
    solids['bone:humerus'] = solids.pop('humerus_GH_ball')
    sources['bone:humerus'] = sources.pop('humerus_GH_ball')
    solids['bone:scapula'] = solids.pop('bone:scapula_GH_only')
    sources['bone:scapula'] = sources.pop('bone:scapula_GH_only')
    intersections = []
    for a, b in itertools.combinations(solids, 2):
        volume = I.intersection_volume(solids[a], solids[b])
        if volume > 1e-9:
            intersections.append(dict(a=a, b=b, volume_mm3=volume,
                category='compliant contact requires deformation review'
                if a.startswith('tpu:') or b.startswith('tpu:') else 'rigid fit requires review'))
    wires = []
    for row in routes:
        wire = D.outer_envelope(np.array(row['full_wire_path_source_mm']), .3)
        hits = [dict(part=k, volume_mm3=I.intersection_volume(wire, s))
                for k, s in solids.items()]
        wires.append(dict(head=row['head'], intersections=hits,
                          hits_over_0_001_mm3=[h for h in hits if h['volume_mm3'] > .001]))
        print('GH-only wire', row['head'], wires[-1]['hits_over_0_001_mm3'], flush=True)
    retained_bone = solids['bone:humerus']
    report = dict(scope='GH-only retained joint comparison, available zero-pose parts; NOT adopted or drive validated',
        n_instances=len(solids), sources=sources,
        input_sha256={p.relative_to(A.ROOT).as_posix(): A.hashed(p)
                      for p in [joint_file, assembly_file, routes_file, Path(__file__).resolve(),
                                Path(A.__file__).resolve()]},
        parent_geometry_sources={'AC': base_meta, 'G1_scapula': eng_meta},
        GH_delta=dict(added_to_AC_context_mm3=float(added.volume()),
                      removed_from_AC_context_mm3=float(removed.volume())),
        anatomical_delta=dict(scapula_removed_mm3=float((oldscap-isolated).volume()),
            scapula_added_mm3=float((isolated-oldscap).volume()),
            humerus_removed_mm3=float((oldhum-retained_bone).volume()),
            humerus_added_mm3=float((retained_bone-oldhum).volume())),
        construction='anatomical scapula - (AC base - G1 scapula) + (G1 scapula - AC base); retain G1 head/cage hardware',
        screw_clearance_reliefs=bolt_reliefs,
        intersections=intersections, wire_probes=wires,
        missing=['real muscle bone-side fasteners, antirotation and crimp fastening',
                 'complete fixture contact and jaws'],
        not_verified=['all motion paths must be recomputed; old G1 passes do not transfer',
            'GH retention under arbitrary escape paths, load and assembly tolerance',
            'TPU contact deformation, real wire force paths and drive range',
            'manufacturing strength and assembly access; screws are unthreaded envelopes'],
        stage_pass=False, actual_drive_range_deg=None)
    (dest/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('GH-only candidate',len(solids),'parts; pairs over .05mm3',
          [r for r in intersections if r['volume_mm3']>.05],flush=True)
    return report


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=DEST)
    run(parser.parse_args().output_dir)
