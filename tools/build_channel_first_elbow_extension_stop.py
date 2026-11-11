"""Removable humeral stop arm and integral forearm lug; explicit zero contact faces."""
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,export,rod
from build_channel_first_inventory import provenance


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'elbow_mount_candidate/assembly.json';parent=read(parent_file)
    definition_file=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/definition.json'
    definition=read(definition_file);E=np.array(definition['centres_source_mm']['elbow']);G=np.array(definition['centres_source_mm']['GH'])
    joint_file=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection/assembly.json'
    joint=read(joint_file);axis=np.array(joint['axis_source']);up=G-E;up-=axis*(up@axis);up/=np.linalg.norm(up)
    side=np.cross(axis,up);F=np.stack([axis,up,side],axis=1)
    def point(v):return E+F@np.array(v)
    def box(lo,hi):return A.I._manifold().Manifold.cube(np.array(hi)-lo).translate(lo).transform(np.c_[F,E])
    def bar(a,b,r):return rod(point(a),point(b),r)
    solids={k:load_component(s) for k,s in parent['sources'].items()}
    hum,fore=solids['bone:humerus'],solids['bone:forearm']
    half=joint['fork_half_gap_mm'];width=joint['ear_width_mm'];outer=half+width
    dest=base/'elbow_extension_stop';dest.mkdir(exist_ok=True)
    trials=[];selected=None
    if (dest/'trial_screen.json').exists():
        import hashlib
        old=(dest/'trial_screen.json').read_bytes()
        (dest/('trial_history_'+hashlib.sha256(old).hexdigest()[:12]+'.json')).write_bytes(old)
    for sign in (1,-1):
        # The fixed arm installs after forearm insertion and pin retention.
        a,b=sorted([sign*outer,sign*(outer+.2)])
        plate_lo,plate_hi=sorted([sign*(outer+.2),sign*(outer+3.2)])
        stock=[];holes=[];back_seats=[];parts={}
        for index,y in enumerate((15.,30.)):
            stock.append(bar([sign*half,y,0],[sign*outer,y,0],4.5))
            holes.append(bar([sign*(half-8),y,0],[sign*(outer+10),y,0],1.7))
            back_seats.append(bar([sign*(half-8),y,0],[sign*half,y,0],3.7))
            washer=bar([sign*outer,y,0],[sign*(outer+.2),y,0],3.5)-holes[-1]
            parts[f'stop:seat_washer_{index}']=washer
            head=outer+3.7
            parts[f'stop:bolt_{index}']=bar([sign*(head-16),y,0],[sign*head,y,0],1.5)+bar([sign*head,y,0],[sign*(head+3),y,0],2.75)
            parts[f'stop:top_washer_{index}']=bar([sign*(outer+3.2),y,0],[sign*head,y,0],3.5)-holes[-1]
            parts[f'stop:back_washer_{index}']=bar([sign*(half-.5),y,0],[sign*half,y,0],3.5)-holes[-1]
            # Hex envelope, explicit through-bore; no thread or preload claim.
            z=axis*sign;x=up;yaxis=np.cross(z,x)
            frame=np.stack([x,yaxis,z],axis=1)
            nut=A.I._manifold().Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[frame,point([sign*(half-2.9),y,0])])
            parts[f'stop:nut_{index}']=nut-holes[-1]
        mf=A.I._manifold();hole_union=mf.Manifold.batch_boolean(holes,mf.OpType.Add)
        pin_void=bar([sign*(half-1),0,0],[sign*(outer+1),0,0],5.85)
        new_hum=(hum+(mf.Manifold.batch_boolean(stock,mf.OpType.Add)-pin_void))-hole_union-mf.Manifold.batch_boolean(back_seats,mf.OpType.Add)
        xlo,xhi=sorted([sign*26.,sign*30.])
        moving_face=box([xlo,-18,6],[xhi,-10,10])
        new_fore=(fore+moving_face+bar([sign*1.5,-7,3],[sign*1.5,-14,8],2.)
                  +bar([sign*1.5,-14,8],[sign*28,-14,8],2.))
        plate=box([plate_lo,-22,0],[plate_hi,36,20])
        for y in (15.,30.):plate+=bar([sign*(outer+.2),y,0],[sign*(outer+3.2),y,0],4.5)
        plate=plate-bar([sign*(outer-1),0,0],[sign*(outer+5),0,0],9.5)
        fixed_face=box([xlo,-18,10],[xhi,-10,14])
        arm=(plate+fixed_face+bar([sign*28,-14,13],[sign*28,-14,18],2.5)+bar([sign*28,-14,18],[sign*(outer+1.7),-14,18],2.5))-hole_union
        # Arm-to-plate bridge ends inside the U plate; connected material is checked.
        parts['stop:arm']=arm
        candidate={**solids,'bone:humerus':new_hum,'bone:forearm':new_fore,**parts}
        changed={'bone:humerus','bone:forearm',*parts};hits=[]
        for ka,kb in combinations(candidate,2):
            if ka not in changed and kb not in changed:continue
            volume=float(A.I.intersection_volume(candidate[ka],candidate[kb]))
            if volume>.001:
                overlap=(candidate[ka]^candidate[kb]).transform(np.c_[F.T,-F.T@E])
                hits.append(dict(left=ka,right=kb,overlap_mm3=volume,
                                 overlap_local_bbox_mm=list(overlap.bounding_box())))
        connected={k:sum(s.volume()>1e-8 for s in candidate[k].decompose()) for k in changed}
        trials.append(dict(side=sign,hits=hits,positive_components=connected))
        if not hits and all(v==1 for v in connected.values()):
            selected=(candidate,parts,sign,moving_face,fixed_face);break
    save(dest/'trial_screen.json',dict(parent=provenance(parent_file),trials=trials))
    if selected is None:
        print('No valid stop candidate; see trial_screen.json',flush=True);return
    candidate,parts,sign,moving_face,fixed_face=selected
    sources=dict(parent['sources'])
    for key in ['bone:humerus','bone:forearm',*parts]:
        _,meta=export(candidate[key].as_original().simplify(1e-6),dest/(key.replace(':','_')+'.obj'),E)
        sources[key]=meta
    save(dest/'assembly.json',dict(sources=sources,display_only_TPU_sources=parent['display_only_TPU_sources'],
        parent=provenance(parent_file),joint_reference=provenance(joint_file),definition=provenance(definition_file),
        generator=provenance(ROOT/'tools/build_channel_first_elbow_extension_stop.py'),frame_columns=F.tolist(),origin_source_mm=E.tolist(),
        selected_side=sign,nominal_contact_angle_deg=0.,stop_face_area_mm2=32.,lever_arm_mm=14.,
        legal_contact=dict(parts=['stop:arm','bone:forearm'],local_plane_z_mm=10.,
            footprint_local_mm=dict(x=sorted([sign*26.,sign*30.]),y=[-18.,-10.]),
            meaning='Zero-gap stop faces only; no whole pair collision exemption'),
        assembly_sequence='Forearm and joint pin installed before removable stop arm; two M3 through-bolts, seat spacers and rear nuts hold arm to added humeral bosses. Full insertion/tools/preload pending.',
        status='zero_pose_stop_candidate_positive_sweep_and_assembly_pending'))
    save(dest/'view_scene.json',dict(sources={**sources,**parent['display_only_TPU_sources']},status='stop_candidate_static_only'))
    print('Stop geometry candidate exported; side',sign,'parts',len(sources),flush=True)


if __name__=='__main__':main()
