"""Readback/contact probes for the fitted proximal mounts, not strength tests."""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import fit_shoulder_origin_mounts as F
import audit_shoulder_retained_joint_fit as A


def run():
    source=F.DEST/'assembly_check.json';assembly=json.loads(source.read_text(encoding='utf8'))
    if assembly['script_sha256']!=A.hashed(Path(F.__file__).resolve()):raise ValueError('Stale fitted assembly')
    solids={};census={}
    for key,row in assembly['sources'].items():
        solids[key],_,meta=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
        census[key]=dict(boundary=meta['boundary_edges'],nonmanifold=meta['nonmanifold_edges'],
            positive_components=sum(c.volume()>1e-6 for c in solids[key].decompose()))
    audit_file=F.B.L.BASE/'origin_fixation_review/fixation_audit.json'
    audit=json.loads(audit_file.read_text(encoding='utf8'));bone=solids['bone:scapula']
    before=json.loads(F.B.SOURCE.read_text(encoding='utf8'))
    rows=[]
    for row in audit['rows']:
        head=row['head'];route_module='infraspinatus' if head=='FJ1500' else 'subscapularis'
        n=np.array(row['holes'][0]['normal_source']);holes=[]
        for hole in row['holes']:
            Q=np.array(hole['axis_point_source_mm'])
            hits=bone.ray_cast(Q+n*10,Q-n*60)
            shaft=F.B.rod(Q+n*10,Q-n*60,1.5)
            holes.append(dict(sign=hole['sign'],remaining_axis_crossings=len(hits),
                M3_shaft_probe_bone_intersection_mm3=A.I.intersection_volume(shaft,bone)))
        cap='cap:'+head+'/origin';tpu='tpu:'+route_module
        front='origin:'+head+'/front_support';back='origin:'+head+'/backing_plate'
        contacts=[]
        for a,b,direction in [(cap,tpu,-n),(front,tpu,n),(front,'bone:scapula',-n),(back,'bone:scapula',n)]:
            contacts.append(dict(a=a,b=b,nominal_material_intersection_mm3=A.I.intersection_volume(solids[a],solids[b]),
                nominal_gap_capped_1mm=A.I.gap(solids[a],solids[b],1.),
                probe_translation_source_mm=(direction*.05).tolist(),
                interference_at_0_05mm_mm3=A.I.intersection_volume(solids[a].translate(direction*.05),solids[b])))
        old_meta=before['sources'][tpu]
        old,_,_=A.load_checked(A.ROOT/old_meta['file'],old_meta['translation_source_mm'],old_meta['sha256'])
        I=np.mean([h['axis_point_source_mm'] for h in row['holes']],axis=0)
        zone=F.B.rod(I-n*8,I+n*.01,11.2);void_changes=[]
        for shell in old.decompose():
            if shell.volume()>=-1e-9:continue
            v,f=A.mesh_geom.manifold_mesh_to_arrays(shell.to_mesh64())
            void=A.I._manifold().Manifold(A.I._manifold().Mesh64(
                np.array(v,dtype=np.float64,order='C',copy=True),
                np.array(f[:,::-1],dtype=np.uint64,order='C',copy=True)))
            added=void^solids[tpu]
            void_changes.append(dict(old_enclosed_void_mm3=float(void.volume()),
                material_now_in_old_void_mm3=float(added.volume()),
                added_outside_proximal_terminal_zone_mm3=float((added-zone).volume()),
                bounds_source_mm=list(void.bounding_box())))
        rows.append(dict(head=head,actual_bone_hole_checks=holes,contact_geometry=contacts,
            closed_void_changes_from_72_object_parent=void_changes,
            interpretation='normal contact geometry only; bolt preload, deformation, friction and load capacity unverified'))
    report=dict(scope=__doc__,assembly_sha256=A.hashed(source),script_sha256=A.hashed(Path(__file__).resolve()),
        generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [Path(F.__file__).resolve(),Path(F.B.__file__).resolve(),
             A.ROOT/'mujoco_arm/tools/fit_shoulder_deltoid_to_collar.py',
             A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py']},
        readback_part_count=len(solids),mesh_census=census,proximal_contact_checks=rows,
        free_terminal_constraints_required=assembly['free_terminal_constraints_required'],
        not_verified=['all real clamp/preload/tool and force-capacity tests','nonlinear flexible body and tendon motion',
                      'existing fixture bone jaws/pads and bench attachment'],stage_pass=False)
    (F.DEST/'mount_contact_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Readback',len(solids),'parts')
    for row in rows:print(row['head'],'holes',row['actual_bone_hole_checks'],
        'contact probes',[(c['a'],c['b'],round(c['nominal_gap_capped_1mm'],5),round(c['interference_at_0_05mm_mm3'],5)) for c in row['contact_geometry']])
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
