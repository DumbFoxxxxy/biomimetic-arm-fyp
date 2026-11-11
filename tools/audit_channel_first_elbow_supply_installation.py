"""Staged inlet assembly with only explicitly not-yet-installed parts absent."""
import numpy as np
import hashlib
import argparse
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--late-flexor-bone0',action='store_true')
    args=parser.parse_args()
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface'
    assembly_file=folder/'assembly.json';assembly=read(assembly_file)
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    stubs={h.removeprefix('sheath:'):load_component(v) for h,v in assembly['conduit_stub_sources'].items()}
    stages=[];rejected_paths=[]
    for row in assembly['interfaces']:
        head=row['id'];axis=np.array(row['reaction_axis_source']);F=np.array(row['frame_columns'])
        prefix='supply:'+head+'_'
        late=[f'elbow:{head}_bone_0_{suffix}' for suffix in ('bolt','front_washer','nut','back_washer')] if args.late_flexor_bone0 and head=='FJ1486' else []
        fixed={k:v for k,v in solids.items() if not k.startswith(prefix) and k not in late}
        def stage(label,moving,vector,distance):
            hits=None;checked=0
            for delta in np.arange(distance,-.001,-.5):
                obstacles=[]
                for key,part in moving.items():
                    placed=part.translate(vector*delta)
                    for target,solid in fixed.items():
                        v=float(A.I.intersection_volume(placed,solid))
                        if v>.001:obstacles.append(dict(moving=key,target=target,overlap_mm3=v))
                checked+=1
                if obstacles:
                    hits=dict(offset_mm=float(delta),interferences=obstacles);break
            stages.append(dict(head=head,stage=label,offset_start_mm=distance,step_mm=.5,
                direction_source=vector.tolist(),checked_states=checked,stationary_keys=list(fixed),
                moving_keys=list(moving),first_interference=hits))
            print(head,label,'checked',checked,'first',hits,flush=True)
            # A failed assembly stage is not followed by fictitious later stages.
            if hits:return False
            fixed.update(moving);return True
        def retainer_docking(key):
            """Search azimuthal docking above the ferrule, not unsupported wire."""
            solid=solids[key];P=np.array(row['origin_source_mm']);tested=0
            fixed_keys=list(fixed);bounds=np.array([fixed[k].bounding_box() for k in fixed_keys])
            def pose_hit(angle,distance):
                placed=solid.transform(A.I.rotation_about(P,axis,angle)).translate(axis*distance)
                b=np.array(placed.bounding_box());admit=np.all(b[:3]<=bounds[:,3:]+1e-8,axis=1)&np.all(bounds[:,:3]<=b[3:]+1e-8,axis=1)
                for j in np.flatnonzero(admit):
                    volume=float(A.I.intersection_volume(placed,fixed[fixed_keys[j]]))
                    if volume>.001:return dict(target=fixed_keys[j],overlap_mm3=volume)
                return None
            for turn in (30.,-30.,60.,-60.,90.,-90.,120.,-120.,180.):
                for height in (2.,4.,6.,8.,10.,12.):
                    tested+=1
                    poses=[(turn,d) for d in np.arange(20.,height-.01,-.5)]
                    poses += [(a,height) for a in np.linspace(turn,0.,round(abs(turn))+1)]
                    poses += [(0.,d) for d in np.arange(height,-.01,-.5)]
                    failed=None
                    for angle,distance in poses:
                        hit=pose_hit(angle,distance)
                        if hit:failed=dict(angle_deg=float(angle),offset_mm=float(distance),hit=hit);break
                    if failed:
                        rejected_paths.append(dict(head=head,kind='retainer_docking',turn_deg=turn,height_mm=height,first=failed));continue
                    # Retainer is prethreaded and held at its parked pose while
                    # ferrule enters. Verify that omitted mutual interference.
                    parked=solid.transform(A.I.rotation_about(P,axis,turn)).translate(axis*20.)
                    for delta in np.arange(20.,-.01,-.5):
                        for item in (solids[prefix+'ferrule'],stubs[head]):
                            if A.I.intersection_volume(parked,item.translate(axis*delta))>.001:
                                failed=dict(parked_plate_vs_entering_ferrule_offset_mm=float(delta));break
                        if failed:break
                    if failed:
                        rejected_paths.append(dict(head=head,kind='parked_retainer',first=failed));continue
                    stages.append(dict(head=head,stage='retainer_azimuthal_docking',
                        parked_azimuth_deg=turn,parking_offset_mm=20.,turnback_height_mm=height,
                        translation_step_mm=.5,rotation_step_deg=1.,checked_states=len(poses),
                        stationary_keys=list(fixed),moving_keys=[key],first_interference=None,
                        alternatives_tested=tested,
                        support='Plate manually held on already threaded long sheath; no driven motion or hidden cable guide.',
                        parked_plate_vs_ferrule_stub_checked=True))
                    fixed[key]=solid;print(head,'retainer docking',turn,height,'clear',flush=True);return True
            return False
        ferrule_key=prefix+'ferrule'
        parked_key='assembly_hold:'+head+'_retainer_at_plus20mm'
        parked=solids[prefix+'retainer'].translate(axis*20.)
        parked_hits=[]
        for target,part in fixed.items():
            volume=float(A.I.intersection_volume(parked,part))
            if volume>.001:parked_hits.append(dict(target=target,overlap_mm3=volume))
        if parked_hits:
            stages.append(dict(head=head,stage='prethreaded_retainer_parking',checked_states=1,
                first_interference=dict(offset_mm=20.,interferences=parked_hits)));continue
        fixed[parked_key]=parked
        if not stage('ferrule_and_nominal_stub_axial_insertion',
                     {ferrule_key:solids[ferrule_key],'nominal_sheath:'+head:stubs[head]},axis,20.):continue
        failed=False
        for i,(a,b) in enumerate(row['ear_centres_local_mm']):
            radial=F@np.array([a,b,0.]);radial/=np.linalg.norm(radial)
            key=prefix+f'nut_{i}'
            if not stage(f'nut_{i}_side_insertion',{key:solids[key]},radial,8.):failed=True;break
        if failed:continue
        key=prefix+'retainer'
        del fixed[parked_key]  # The same held plate now moves; no duplicate part.
        if not stage('retainer_prethreaded_on_sheath_then_axial_seating',{key:solids[key]},axis,20.):
            rejected_paths.append(stages.pop())
            if not retainer_docking(key):continue
        for i in (0,1):
            key=prefix+f'head_washer_{i}'
            if not stage(f'head_washer_{i}_axial_seating',{key:solids[key]},axis,5.):failed=True;break
        if failed:continue
        for i in (0,1):
            key=prefix+f'bolt_{i}'
            if not stage(f'M2_{i}_axial_insertion_envelope',{key:solids[key]},axis,15.):failed=True;break
        if failed:continue
        if late:
            mounts=read(ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_mount_candidate/assembly.json')
            anchor=next(r for r in mounts['mount_seats'] if r['muscle']==head and r['anchor']==0)
            normal=np.array(anchor['axis_source'])
            front={k:solids[k] for k in late if k.endswith(('bolt','front_washer'))}
            rear={k:solids[k] for k in late if k.endswith(('nut','back_washer'))}
            if not stage('late_bone0_front_bolt_and_washer',front,normal,45.):continue
            stage('late_bone0_rear_nut_and_washer',rear,-normal,35.)
    output=folder/('installation_late_bone0.json' if args.late_flexor_bone0 else 'installation_screen.json')
    if output.exists():
        old=output.read_bytes()
        (folder/('installation_history_'+hashlib.sha256(old).hexdigest()[:12]+'.json')).write_bytes(old)
    save(output,dict(assembly=provenance(assembly_file),stages=stages,rejected_paths=rejected_paths,
        late_flexor_bone0=args.late_flexor_bone0,
        generator=provenance(ROOT/'tools/audit_channel_first_elbow_supply_installation.py'),
        all_planned_stages_completed=len(stages)==(18 if args.late_flexor_bone0 else 16) and all(s['first_interference'] is None for s in stages),
        scope='Zero-pose subassembly screening. Original cap/bone frame at source pose; other inlet assembled. Optional flexor bone0 hardware explicitly absent until late front/rear insertion, so frame needs assembly holding, not driven loading. Retainer prethreaded before ferrule termination. Steel core not yet threaded. Short sheath reference only. No whole-neighbour exemptions. Threads are envelopes; tools, preload, manual holding design, crimp/bond, original muscle insertion and full supply hose absent.'))


if __name__=='__main__':main()
