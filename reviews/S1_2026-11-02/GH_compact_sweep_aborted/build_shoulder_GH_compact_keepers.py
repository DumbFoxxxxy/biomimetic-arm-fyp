"""Rebuilt recessed GH keepers with a local flexion-motion relief.

Starts at the 187-object entry candidate. Recesses heads/nuts, moves keeper 0
inward 0.3 mm, and trims only added pocket material around sampled bare-bone
flexion. This is nominal geometry, not load/strength or muscle-drive proof.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_bench_mount as B
from open_shoulder_GH_assembly_entries import DEST as PARENT
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_GH_assembly_review/compact_recessed_keepers'


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    source=PARENT/'assembly_check.json';parent=json.loads(source.read_text(encoding='utf8'));sources=dict(parent['sources'])
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in sources.items()}
    reference=A.G1/'assembly.json';joint=json.loads(reference.read_text(encoding='utf8'));G=np.array(joint['centre_source_mm']);n=np.array(joint['neck_axis']);mf=A.I._manifold()
    old=S['bone:scapula'];body=old;drills=[];definitions=[]
    for i,entry in enumerate(parent['GH_assembly_entries']['definitions'][:2]):
        u=np.array(entry['direction_source']);v=np.cross(n,u);F=np.column_stack([u,v,n]);c=np.array(joint['bolt_sites_display_mm'][i])+A.frames.SHIFT
        station=4.4 if i==0 else 4.7
        boss=B.box([-4.6,-5.5,-8.5],[7.2,5.5,-3.])
        boss-=B.box([-3.35,-2.85,-7.3],[8.,2.85,-4.7]);boss-=B.box([-3.65,-3.6,-4.8],[8.,3.6,-3.95])
        boss-=rod([0,0,-10],[0,0,1],1.7)
        boss=boss.transform(np.c_[F,c])-mf.Manifold.sphere(joint['socket_inner_radius_mm'],40).translate(G)
        body+=boss
        pin=c+F@np.array([station,0.,-5.2]);w=A.D.unit(np.cross(np.eye(3)[np.argmin(abs(v))],v));H=np.column_stack([w,np.cross(v,w),v])
        drill=rod(pin-v*12.,pin+v*12.,1.1)+rod(pin-v*25.,pin-v*3.5,2.05)
        drill+=mf.Manifold.cylinder(21.6,4.2/np.sqrt(3),4.2/np.sqrt(3),6).transform(np.c_[H,pin+v*3.4]);drills.append(drill)
        screw=rod(pin-v*5.5,pin-v*3.5,1.9)+rod(pin-v*3.5,pin+v*5.5,1.)
        screw-=mf.Manifold.cylinder(1.2,1.55/np.sqrt(3),1.55/np.sqrt(3),6).transform(np.c_[H,pin-v*5.6])
        nut=mf.Manifold.cylinder(1.6,4./np.sqrt(3),4./np.sqrt(3),6)-mf.Manifold.cylinder(1.8,1.1,1.1,48).translate([0,0,-.1])
        nut=nut.transform(np.c_[H,pin+v*3.5]);hardware={}
        for name,part in [('screw',screw),('nut',nut)]:
            key=f'GH_keeper:{i}/{name}';actual,meta=local_export(part,DEST/f'compact_keeper_{i}_{name}.stl',c)
            if not meta['mesh_checks_pass']:raise ValueError('Invalid compact fastener '+key)
            S[key]=actual;sources[key]=meta;hardware[name]=key
        definitions.append(dict(index=i,source_origin_mm=c.tolist(),basis_columns_source=F.tolist(),entry_direction_source=u.tolist(),
                                keeper_axis_source=v.tolist(),cross_pin_local_u_n_mm=[station,-5.2],
                                keeper_hardware=hardware,keyed_nut_channel_width_mm=5.7,nominal_M3_nut_across_flats_mm=5.5,
                                washer_channel_width_mm=7.2,cross_bore_diameter_mm=2.2,
                                candidate_screw_under_head_length_mm=9.,candidate_M2_head_height_mm=2.,candidate_head_radius_mm=1.9,
                                head_recess_bearing_local_v_mm=-3.5,head_top_local_v_mm=-5.5,
                                nut_recess_across_flats_mm=4.2,nut_local_v_interval_mm=[3.5,5.1],
                                socket_across_flats_mm=1.55,socket_depth_mm=1.1,
                                nominal_axis_shift_from_external_keeper_mm=-.3 if i==0 else 0.))
    for drill in drills:body-=drill
    # Bound the design keepout to the local GH region; original entry-bone
    # material is never removed by this motion-profile operation.
    crop=mf.Manifold.sphere(40.,48).translate(G)
    root=S['bone:humerus']^crop;axis=A.frames.gh_probe_axes()['flexion']
    angles=np.linspace(80.,120.,81)
    posed=[root.transform(A.I.rotation_about(G,axis,float(a))) for a in angles]
    envelope=mf.Manifold.batch_boolean(posed,mf.OpType.Add).minkowski_sum(mf.Manifold.cube([.4]*3,True))
    new_material=body-old;trim=new_material^envelope;body-=trim
    print('Added pocket material relieved',trim.volume(),flush=True)
    if sum(p.volume()>1e-6 for p in body.decompose())!=1:raise ValueError('Compact relief disconnects engineered scapula')
    origin=np.array(sources['bone:scapula']['translation_source_mm']);actual,meta=local_export(body,DEST/'scapula_compact_keepers.stl',origin)
    if not meta['mesh_checks_pass']:raise ValueError('Compact scapula readback fails')
    S['bone:scapula']=actual;sources['bone:scapula']=meta
    overlaps=[]
    for a,b in itertools.combinations(S,2):
        volume=A.I.intersection_volume(S[a],S[b])
        if volume>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=volume,compliant=is_compliant(a) or is_compliant(b)))
    routes=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';wires=[]
    for r in json.loads(routes.read_text(encoding='utf8'))['routes']:
        wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3);hits={k:A.I.intersection_volume(wire,s) for k,s in S.items()}
        wires.append(dict(head=r['head'],hits_over_threshold=[dict(part=k,volume_mm3=z) for k,z in hits.items() if z>.001]))
    report=dict(scope=__doc__,n_instances=len(S),sources=sources,parent_report=source.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),joint_reference_sha256=A.hashed(reference),routes_sha256=A.hashed(routes),
                generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                    [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py',A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py']},
                GH_nut_keepers=dict(definitions=definitions,scapula_changes=dict(added_mm3=float((actual-old).volume()),removed_mm3=float((old-actual).volume())),
                                   motion_relief=dict(added_material_removed_mm3=float(trim.volume()),source_angles_deg=angles.tolist(),
                                                     source_crop_radius_mm=40.,xyz_each_side_keepout_mm=.2,
                                                     original_bone_material_intentionally_removed_by_motion_relief_mm3=0.)),
                GH_assembly_entries=parent['GH_assembly_entries'],GH_retention_revision=parent['GH_retention_revision'],GH_lip_relief=parent['GH_lip_relief'],
                free_terminal_constraints_required=parent['free_terminal_constraints_required'],
                intended_body_owners={**parent['intended_body_owners'],**{k:'fixed_context' for k in sources if k.startswith('GH_keeper:')}},
                pair_count_checked=len(S)*(len(S)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['bare joint motion and staged assembly after readback','reaction faces, torque, seat thickness and physical retention',
                              'candidate screw selection/length, thread, actual tools and preload','full muscle assembly, deformable motion and drive'])
    (DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Changes',report['GH_nut_keepers']['scapula_changes'],flush=True)
    print('Rigid overlaps',[r for r in overlaps if r['volume_mm3']>.001 and not r['compliant']],flush=True)
    print('Wire failures',[r for r in wires if r['hits_over_threshold']],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
