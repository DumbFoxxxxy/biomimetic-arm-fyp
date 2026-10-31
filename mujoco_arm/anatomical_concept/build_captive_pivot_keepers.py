"""Two thin U keepers capture four grooved pin heads with real cross-fasteners.

Metal keeper candidate; thin sheet forming and fastener procurement dimensions
remain unqualified. No blind threading or friction-only pin fixation is assumed.
"""
from build_direct_muscle_yoke import *
from build_forearm_receivers import bolt
from pathlib import Path

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';dest=parent/'captive_keepers';dest.mkdir(exist_ok=True);a=json.loads((parent/'assembly.json').read_text());rows={r['key']:dict(r) for r in a['instances']};original={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};motion=json.loads((parent/'conditional_motion.json').read_text())['poses'];module={};sources=[parent/'assembly.json',parent/'conditional_motion.json']+[ROOT/r['file'] for r in rows.values()]
    for fid,folder in [('FJ1486','cap_swivel_flexor'),('FJ1477','cap_swivel_study')]:
        p=DEST/folder/'screen.json';s=json.loads(p.read_text());C=np.array(s['pivot_zero_display_mm']);e=np.array(s['trunnion_axis_zero_display']);choice=s['selected_fork'];z=np.array(json.loads((muscle_folder(fid)/'design.json').read_text())['basis_columns_display'])[:,2];delta=np.array(choice.get('bridge_delta_display_mm',np.array(choice['bridge_direction_display'])*choice.get('bridge_reach_mm',20)));direction=unit(delta);L=np.linalg.norm(delta);normal=unit(np.cross(direction,e));F=np.column_stack([e,normal,direction]);module[fid]=(C,e,direction,L,F);sources += [p,muscle_folder(fid)/'design.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};refs={}
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in rows:
            p=OUT/'repaired_supports'/(r['id']+'.obj');refs['reference_'+r['id']]=load(p,raw=True)[0];hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    trials=[];selected=None
    for off_flex,off_ext in itertools.product([8,-8,0],[8,-8,0]):
        frame=original['paired_bone_frame'];parts={};features=[];bosses=[];clamp_holes=[]
        for fid,offset in [('FJ1486',off_flex),('FJ1477',off_ext)]:
            C,e,direction,L,F=module[fid];T=np.c_[F,C];keeper_pieces=[]
            for sign in [-1,1]:
                lo,hi=sorted([sign*19.15,sign*19.55]);plate=rod(C+e*lo,C+e*hi,2.5,48)+mf.Manifold.cube([.4,5,L+3.4]).translate([lo,-2.5,0]).transform(T)
                plate=plate-rod(C+e*(lo-.2),C+e*(hi+.2),1.65,48)-mf.Manifold.cube([.8,3.3,5]).translate([lo-.2,-1.65,-5]).transform(T);keeper_pieces.append(plate)
                key=fid+'_pin_'+str(sign);pin=original[key];v=e*sign;groove=rod(C+v*19.05,C+v*19.65,3,48)-rod(C+v*18.9,C+v*19.8,1.6,48);parts[key]=(pin-groove).as_original().simplify(.0001)
            keeper_pieces.append(mf.Manifold.cube([39.1,5,.4]).translate([-19.55,-2.5,L+3.0]).transform(T));keeper=mf.Manifold.batch_boolean(keeper_pieces,mf.OpType.Add)
            q=C+e*offset+direction*L;boss=mf.Manifold.cube([6,5,5.5]).translate([offset-3,-2.5,L-2.5]).transform(T);hole=rod(q-direction*6,q+direction*6,1.1,48);bosses.append(boss);clamp_holes.append(hole);frame=(frame+boss)-hole;parts[fid+'_keeper']=(keeper-hole).as_original().simplify(.001)
            fastener_frame=F;Tf=np.c_[fastener_frame,q]
            parts[fid+'_keeper_bolt']=bolt(10,3.9,1,2,1.3,32).transform(Tf)
            parts[fid+'_keeper_front_washer']=(rod(q+direction*3.4,q+direction*3.9,2.5)-hole)
            parts[fid+'_keeper_back_washer']=(rod(q-direction*2.8,q-direction*2.5,2.5)-hole)
            nut=mf.Manifold.cylinder(1.6,4/np.sqrt(3),4/np.sqrt(3),6).translate([0,0,-4.4]).transform(Tf)-hole;parts[fid+'_keeper_nut']=nut
            features.append({'muscle':fid,'keeper_bridge_direction_display':direction.tolist(),'keeper_bolt_centre_display_mm':q.tolist(),'keeper_bolt_axis_display':direction.tolist(),'keeper_bolt_trunnion_offset_mm':offset,'nominal_clamp_bolt_length_mm':10,'sheet_thickness_mm':.4,'pin_head_groove_axial_interval_mm':[19.05,19.65],'groove_root_radius_mm':1.6,'keeper_inner_radius_mm':1.65,'keeper_axial_interval_mm':[19.15,19.55]})
        unrelieved=frame
        for fid in module:
            direction=module[fid][2];sweep=parts[fid+'_keeper'].minkowski_sum(rod([0,0,0],direction*5,.01,8));frame=frame-sweep
        # Re-establish the exact planar clamp faces and bores. The keeper
        # moves away from these faces during insertion, so they do not obstruct
        # the path. Do not replace a positive bearing face with a hidden gap.
        frame=(frame+mf.Manifold.batch_boolean(bosses,mf.OpType.Add))-mf.Manifold.batch_boolean(clamp_holes,mf.OpType.Add)
        relieved=max(0,float((unrelieved-frame).volume()));assert relieved<100,'plate relief exceeds local correction scope'
        frame=frame.as_original().simplify(.0001);parts['paired_bone_frame']=frame;world={**original,**parts};hard=[];soft=[];components={k:len([s for s in v.decompose() if s.volume()>.001]) for k,v in parts.items()}
        for key,s in parts.items():
            for other,body in world.items():
                if key==other or other in parts and other<key:continue
                v=intersection(s,body)
                if v>.001:(soft if rows.get(other,{}).get('kind')=='TPU' else hard).append({'pair':[key,other],'volume_mm3':v})
        first=None
        if not hard and not soft and all(n==1 for n in components.values()):
            for pose in motion:
                Tf=np.array(pose['forearm_pose']);hits=[]
                for key,s in parts.items():
                    current=s.transform(Tf)
                    for other,body in world.items():
                        if other==key or other in parts or rows.get(other,{}).get('kind')=='TPU':continue
                        owner=rows[other]['owner'];Tother=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;target=body.transform(Tother) if Tother is not None else body;v=intersection(current,target)
                        if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
                    for other,body in refs.items():
                        v=intersection(current,body)
                        if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
                if hits:first={'elbow_deg':pose['elbow_deg'],'hits':hits};break
        trial={'flexor_offset_mm':off_flex,'extensor_offset_mm':off_ext,'local_frame_plate_relief_mm3':relieved,'neutral_rigid_interferences':hard,'neutral_TPU_review_flags':soft,'positive_components':components,'first_conditional_interference':first};trials.append(trial);print('keepers',off_flex,off_ext,'hard',len(hard),'soft',len(soft),'first',first,flush=True)
        if not hard and not soft and first is None and all(n==1 for n in components.values()):selected=(parts,features,trial);break
    report={'scope':'real grooved-head capture and through-bolted Ukeepers, changed geometry vs full nominal parent arm, conditional prescribed cap angles; strength/printing/passive equilibrium excluded','source_sha256':hashes,'trials':trials,'selected':selected[2] if selected else None,'physical_retention_or_looseness_verified':False,'installed_assembly_paths_verified':False,'frozen':False}
    if selected:
        parts,features,trial=selected;capture=[];metas=[]
        for fid,(C,e,direction,L,F) in module.items():
            for sign in [-1,1]:
                key=fid+'_pin_'+str(sign);pin=parts[key];keeper=parts[fid+'_keeper'];v=e*sign;capture.append({'pin':key,'rest_keeper_overlap_mm3':intersection(pin,keeper),'outward_0_09mm_keeper_overlap_mm3':intersection(pin.translate(v*.09),keeper),'outward_0_11mm_keeper_overlap_mm3':intersection(pin.translate(v*.11),keeper),'inward_0_01mm_washer_overlap_mm3':intersection(pin.translate(-v*.01),world[fid+'_washer_'+str(sign)]),'minimum_stub_engagement_at_outward_stop_mm':3.2})
        report.update(features=features,capture_witnesses=capture,conditional_samples=len(motion))
        for key,s in parts.items():
            metas.append(write(dest/(key+'.stl'),s));rows[key]={'key':key,'file':str((dest/(key+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':'elbow'}
        out=dict(a);out.update(scope='paired swivel with actual Ukeeper pin capture, fabrication andloading unqualified',instances=list(rows.values()),frozen=False);(dest/'assembly.json').write_text(json.dumps(out,indent=2),encoding='utf8');(dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':metas},indent=2),encoding='utf8');print('keeper candidate',len(rows),'instances',len(metas),'new/revisedCADparts; capture',capture,flush=True)
    (dest/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())

if __name__=='__main__':run()
