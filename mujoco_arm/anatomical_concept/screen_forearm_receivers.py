"""M8 installed receiver direction probes, not a completed tendon route.

Upstream references are M7 sheave centres, NOT steel contact/control points.
No cable bend or constraint is imposed there. This deliberately diagnoses the
old lever anchor phases before inventing unsupported downstream supports.
"""
from build_forearm_receivers import *
import sys
SCREEN_DIR=DEST/'rephased_anchors' if '--rephase' in sys.argv else DEST

def phase_lever(definition,targets):
    E=np.array(definition['joint_centre_display_mm']);axis=np.array(definition['joint_axis']);n=-axis;down=np.array(definition['zero_forearm_direction']);side=np.cross(axis,down);centre=E+n*56;spine=centre+down*15;endspine=centre+down*58;points={};phases={};ears=[];roots=[]
    for fid,label in [('FJ1486','flexor'),('FJ1477','extensor')]:
        delta=targets[fid]-E;theta=np.degrees(np.arctan2(delta@side,delta@down));phase=(theta-165 if label=='flexor' else theta+15)%360;phases[label]=float(phase);points[label]=centre+15*(down*np.cos(np.radians(phase))+side*np.sin(np.radians(phase)))
    for anchor in definition['bone_anchors']:
        face=np.array(anchor['front_spotface_mm']);root=face+n*1.65;roots.append(root);ears.append(rod(face+n*.15,face+n*3.15,5)-rod(face-n*8,face+n*8,1.7,32))
    end=points['flexor'];opposite=points['extensor'];pieces=ears+[rod(spine,endspine,3),rod(opposite,spine,3),rod(end,spine,3)]
    for root in roots:
        turn=root+down*8;outer=turn+n*((centre-turn)@n);path=[root,turn,outer,end];pieces.extend(rod(a,b,3) for a,b in zip(path,path[1:]))
    lever=mf.Manifold.batch_boolean(pieces,mf.OpType.Add)
    for P in points.values():lever=lever+rod(P-axis*2,P+axis*2,5)
    for anchor in definition['bone_anchors']:
        p=np.array(anchor['front_spotface_mm']);lever=lever-rod(p-n*2,p+n*20,1.7,32)-rod(p+n*3,p+n*20,3.75,32)
    return lever,points,phases

def receiver_pose(parts,standoff,F,P,target):
    local=F.T@(target-P);yaw=np.degrees(np.arctan2(local[1],local[0]));heading=np.array([np.cos(np.radians(yaw)),np.sin(np.radians(yaw)),0]);B=np.array([0.,0.,5.5+standoff]);H=B+8*heading;delta=local-H;pitch=np.degrees(np.arctan2(delta[2],np.linalg.norm(delta[:2])));moving=articulate(parts,yaw,pitch,standoff);world={k:s.transform(np.c_[F,P]) for k,s in moving.items()};direction=heading*np.cos(np.radians(pitch))+np.array([0,0,1])*np.sin(np.radians(pitch));outlet=P+F@(H+direction*9.8)
    return world,float(yaw),float(pitch),outlet

def run():
    a=json.loads((PARENT/'assembly.json').read_text());rows={r['key']:r for r in a['instances']}
    support_path=ROOT/'骨骼/支撑参考/sources.json';support=json.loads(support_path.read_text(encoding='utf8'));references_to_add={r['id']:OUT/'repaired_supports'/(r['id']+'.obj') for r in support['models'] if r['kind']=='bone_reference' and r['id'] not in rows}
    source_paths=[PARENT/'assembly.json',BASE/'design.json',OUT/'elbow_guides/paired_placement/selection.json',support_path]+[ROOT/r['file'] for r in rows.values()]+list(references_to_add.values())
    for fid in ['FJ1486','FJ1477']:
        source_paths.append(DEST/fid/'design.json');source_paths.extend(DEST/fid/(k+'.obj') for k in json.loads((DEST/fid/'design.json').read_text())['parts'])
    source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    original={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};definition=json.loads((BASE/'design.json').read_text());E=np.array(definition['joint_centre_display_mm']);axis=np.array(definition['joint_axis']);n=-axis;down=unit(np.array(definition['zero_forearm_direction']));F0=np.stack([down,np.cross(n,down),n],axis=1)
    for fid,p in references_to_add.items():
        key='reference_'+fid;original[key]=load(p,raw=True)[0];rows[key]={'key':key,'owner':'world','kind':'bone_reference'}
    SCREEN_DIR.mkdir(exist_ok=True)
    references=json.loads((OUT/'elbow_guides/paired_placement/selection.json').read_text())['placements'];lever=original['side_lever'];before=lever;points={k:np.array(v) for k,v in definition['eyelets_display_mm'].items()};phases=None
    if '--rephase' in sys.argv:lever,points,phases=phase_lever(definition,{fid:np.array(v['centre_display_mm']) for fid,v in references.items()})
    for label in ['flexor','extensor']:
        P=points[label];lever=lever-rod(P-n*10,P+n*10,1.7,32)-rod(P+n*2,P+n*20,3.9,32)-rod(P-n*20,P-n*2,3.65,32)
    lever=lever.as_original().simplify(.001);row=write(SCREEN_DIR/'lever_receiving_spotfaces.stl',lever);assert row['components']==1;(SCREEN_DIR/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[row]}),encoding='utf8')
    data={};world_base={k:s for k,s in original.items() if k not in ['side_lever','flexor_tendon_stop','extensor_tendon_stop']};scope=set(json.loads((PARENT/'rigid_motion_scope.json').read_text())['included'])|{'reference_'+fid for fid in references_to_add};neutral=[];pose_rows=[];first=None
    for fid,label in [('FJ1486','flexor'),('FJ1477','extensor')]:
        d=json.loads((DEST/fid/'design.json').read_text());data[fid]=(d,{k:load(DEST/fid/(k+'.obj'))[0] for k in d['parts']},points[label],np.array(references[fid]['centre_display_mm']))
    scalar_cache={}
    def scalar_pose(angle):
        if angle not in scalar_cache:
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();values=[angle]
            for fid,(d,parts,P0,target) in data.items():
                _,yaw,pitch,_=receiver_pose({},d['standoff_mm'],R@F0,E+R@(P0-E),target);values.extend([yaw,pitch])
            scalar_cache[angle]=np.array(values)
        return scalar_cache[angle]
    def angular_step(a,b):
        delta=scalar_pose(b)-scalar_pose(a);delta[[1,3]]=(delta[[1,3]]+180)%360-180;return float(abs(delta).max())
    angles=[0.]
    def add_interval(a,b):
        if angular_step(a,b)>1.00000001:
            mid=(a+b)/2;add_interval(a,mid);add_interval(mid,b)
        else:angles.append(float(b))
    for angle in range(150):add_interval(float(angle),float(angle+1))
    for angle in angles:
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];others={k:(s.transform(T) if rows[k]['owner']=='elbow' else s) for k,s in world_base.items() if k in scope};others['revised_lever']=lever.transform(T);all_receivers={};metrics={}
        for fid,(d,parts,P0,target) in data.items():
            P=E+R@(P0-E);F=R@F0;posed,yaw,pitch,outlet=receiver_pose(parts,d['standoff_mm'],F,P,target)
            for key,s in posed.items():all_receivers[fid+'_'+key]=s
            arm=outlet-E;force=unit(target-outlet);moment=float(np.cross(arm,force)@axis)
            metrics[fid]={'yaw_relative_lever_deg':yaw,'pitch_deg':pitch,'reference_line_moment_arm_mm':moment,'terminal_outlet_display_mm':outlet.tolist(),'reference_point_display_mm':target.tolist(),'pitch_exceeds_local_sampled_80deg':abs(pitch)>80}
        h=[]
        for other,b in others.items():
            if other=='revised_lever':continue
            v=intersection(others['revised_lever'],b)
            if v>.001:h.append({'pair':['revised_lever',other],'volume_mm3':v})
        for key,s in all_receivers.items():
            for other,b in others.items():
                v=intersection(s,b)
                if v>.001:h.append({'pair':[key,other],'volume_mm3':v})
        for (ka,sa),(kb,sb) in itertools.combinations(all_receivers.items(),2):
            if ka[:6]==kb[:6]:continue
            v=intersection(sa,sb)
            if v>.001:h.append({'pair':[ka,kb],'volume_mm3':v})
        if h and first is None:first={'angle_deg':angle,'hits':h}
        pose_rows.append({'elbow_angle_deg':angle,'modules':metrics,'rigid_interferences':h})
        if angle==0:
            for other,b in world_base.items():
                v=intersection(lever,b)
                if v>.001:neutral.append({'pair':['revised_lever',other],'volume_mm3':v,'kind':'TPU_review' if rows[other]['kind']=='TPU' else 'rigid'})
            for key,s in all_receivers.items():
                write(SCREEN_DIR/(key+'_neutral.stl'),s)
                for other,b in world_base.items():
                    v=intersection(s,b)
                    if v>.001:neutral.append({'pair':[key,other],'volume_mm3':v,'kind':'TPU_review' if rows[other]['kind']=='TPU' else 'rigid'})
    summary=[]
    for fid in data:
        module_poses=[r['modules'][fid] for r in pose_rows];moment=[r['reference_line_moment_arm_mm'] for r in module_poses];summary.append({'muscle':fid,'pitch_min_deg':min(r['pitch_deg'] for r in module_poses),'pitch_max_deg':max(r['pitch_deg'] for r in module_poses),'reference_moment_min_mm':min(moment),'reference_moment_max_mm':max(moment),'reference_flexion_sign_count':sum(m>0 for m in moment),'reference_extension_sign_count':sum(m<0 for m in moment)})
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in source_hashes.items()),'Input geometry changed during screen'
    report={'scope':'prescribed receiving-gimbal direction poses referenced to existing sheave CENTRES, not steel contact points or real downstream supports; no complete cable constraint, motor/TPU stroke or drive proof','samples':len(angles),'maximum_elbow_step_deg':float(np.diff(angles).max()),'maximum_step_including_receiver_hinges_deg':max(angular_step(a,b) for a,b in zip(angles,angles[1:])),'return_same_prescribed_rigid_configurations':True,'first_rigid_limit':first,'neutral_contacts':neutral,'modules':summary,'poses':pose_rows,'candidate_phases_deg':phases,'candidate_eyelets_display_mm':{k:p.tolist() for k,p in points.items()},'lever_source_sha256':hashlib.sha256((ROOT/rows['side_lever']['file']).read_bytes()).hexdigest(),'lever_material_removed_mm3':float((before-lever).volume()),'lever_material_added_mm3':float((lever-before).volume()),'source_sha256':source_hashes,'old_uncrimped_stops_removed_from_this_candidate_only':True,'full_transmission_exists':False,'frozen':False};(SCREEN_DIR/'mounted_direction_screen.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('M8 lever components',row['components'],'samples',len(angles),'first rigid limit',first,'neutral flags',len(neutral),flush=True);print('direction probes',summary,flush=True)

if __name__=='__main__':run()
