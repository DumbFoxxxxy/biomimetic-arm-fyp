"""Recess GH cross-keeper screw heads and keyed nuts inside existing cheeks.

Retains M3 keyed pockets and pin centre lines. Removes four keeper washers,
uses recessed candidate M2 cap screws and hex nut seats. No selected catalogue
fastener, material, preload, thread engagement strength or drive acceptance.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_GH_nut_keepers import DEST as PARENT
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_GH_assembly_review/recessed_keeper_candidate'


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    source=PARENT/'assembly_check.json';parent=json.loads(source.read_text(encoding='utf8'))
    sources={k:r for k,r in parent['sources'].items() if not(k.startswith('GH_keeper:') and k.endswith(('/near_washer','/far_washer')))}
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in sources.items()}
    body=S['bone:scapula'];old=body;mf=A.I._manifold();definitions=[]
    for d in parent['GH_nut_keepers']['definitions']:
        i=d['index'];c=np.array(d['source_origin_mm']);F=np.array(d['basis_columns_source']);u,v,n=F.T
        pin=c+F@np.array([4.7,0.,-5.2])
        w=A.D.unit(np.cross(np.eye(3)[np.argmin(abs(v))],v));H=np.column_stack([w,np.cross(v,w),v])
        # Cap-head bearing plane at local v=-3.5; top is flush at -5.5.
        body-=rod(pin-v*5.6,pin-v*3.5,2.05)
        # The hex opening continues through the outer face for axial loading.
        hex_bore=mf.Manifold.cylinder(21.6,4.2/np.sqrt(3),4.2/np.sqrt(3),6).transform(np.c_[H,pin+v*3.4])
        body-=hex_bore
        screw=rod(pin-v*5.5,pin-v*3.5,1.9)+rod(pin-v*3.5,pin+v*5.5,1.)
        socket=mf.Manifold.cylinder(1.2,1.55/np.sqrt(3),1.55/np.sqrt(3),6).transform(np.c_[H,pin-v*5.6])
        screw-=socket
        nut=mf.Manifold.cylinder(1.6,4./np.sqrt(3),4./np.sqrt(3),6)-mf.Manifold.cylinder(1.8,1.1,1.1,48).translate([0,0,-.1])
        nut=nut.transform(np.c_[H,pin+v*3.5])
        hardware={}
        for name,part in [('screw',screw),('nut',nut)]:
            key=f'GH_keeper:{i}/{name}';actual,meta=local_export(part,DEST/f'recessed_keeper_{i}_{name}.stl',c)
            if not meta['mesh_checks_pass']:raise ValueError('Invalid recessed fastener '+key)
            S[key]=actual;sources[key]=meta;hardware[name]=key
        definitions.append({**d,'keeper_hardware':hardware,'candidate_screw_under_head_length_mm':9.,
                            'nominal_M2_screw_length_mm':9.,'candidate_M2_head_height_mm':2.,
                            'candidate_screw_head_height_mm':2.,'candidate_head_top_local_v_mm':-5.5,
                            'screw_shank_tip_local_v_mm':5.5,'nut_local_v_interval_mm':[3.5,5.1],
                            'head_recess_radius_mm':2.05,'head_recess_bearing_local_v_mm':-3.5,
                            'nut_recess_across_flats_mm':4.2,'socket_across_flats_mm':1.55,
                            'socket_depth_mm':1.1,'previous_washer_keys_removed':[f'GH_keeper:{i}/near_washer',f'GH_keeper:{i}/far_washer']})
    if sum(p.volume()>1e-6 for p in body.decompose())!=1:raise ValueError('Recesses disconnect engineered scapula')
    origin=np.array(sources['bone:scapula']['translation_source_mm']);actual,meta=local_export(body,DEST/'scapula_recessed_keepers.stl',origin)
    if not meta['mesh_checks_pass']:raise ValueError('Recessed scapula readback fails')
    S['bone:scapula']=actual;sources['bone:scapula']=meta
    overlaps=[]
    for a,b in itertools.combinations(S,2):
        volume=A.I.intersection_volume(S[a],S[b])
        if volume>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=volume,compliant=is_compliant(a) or is_compliant(b)))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';wires=[]
    for r in json.loads(routefile.read_text(encoding='utf8'))['routes']:
        wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)
        hits={k:A.I.intersection_volume(wire,s) for k,s in S.items()}
        wires.append(dict(head=r['head'],hits_over_threshold=[dict(part=k,volume_mm3=z) for k,z in hits.items() if z>.001]))
    report=dict(scope=__doc__,n_instances=len(S),sources=sources,parent_report=source.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),routes_sha256=A.hashed(routefile),
                GH_nut_keepers=dict(definitions=definitions,scapula_changes=dict(removed_mm3=float((old-actual).volume()),added_by_readback_mm3=float((actual-old).volume()))),
                GH_assembly_entries=parent['GH_assembly_entries'],GH_retention_revision=parent['GH_retention_revision'],GH_lip_relief=parent['GH_lip_relief'],
                free_terminal_constraints_required=parent['free_terminal_constraints_required'],
                intended_body_owners={k:r for k,r in parent['intended_body_owners'].items() if k in sources},
                pair_count_checked=len(S)*(len(S)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
                removed_objects=[k for k in parent['sources'] if k not in sources],
                generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                    [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py']},
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['bare joint motion and staged assembly must be independently rechecked',
                              'candidate screw length/catalogue availability, real thread, tightening and tools',
                              'thin recessed seat strength, print tolerances, wear and physical holding',
                              'full muscle assembly, coupled flexible motion and actuation'])
    (DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Scapula changes',report['GH_nut_keepers']['scapula_changes'],flush=True)
    print('Rigid overlaps',[r for r in overlaps if r['volume_mm3']>.001 and not r['compliant']],flush=True)
    print('Wire failures',[r for r in wires if r['hits_over_threshold']],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
