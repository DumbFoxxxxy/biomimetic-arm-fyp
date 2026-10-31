"""Real bone-rooted common support for two passive distal cap pivots.

Angles are prescribed conditional geometry, not solved passive equilibrium.
The old integrated yoke is replaced by a connected frame on the same seats;
it is not omitted from a purported installed clearance test.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation
from project_runtime import WorkspacePath as Path

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST;dest=parent/'paired_swivel_support';dest.mkdir(exist_ok=True);a=json.loads((parent/'assembly.json').read_text());rows={r['key']:dict(r) for r in a['instances'] if r['key']!='direct_muscle_yoke'};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);down=unit(np.array(j['zero_forearm_direction']));n=np.array(j['bone_anchors'][0]['bolt_axis']);side=unit(np.cross(n,down));modules={};newmeta=[];fixed_new={};forks=[];voids=[];sources=[parent/'assembly.json',BASE/'design.json']+[ROOT/r['file'] for r in rows.values()]
    for fid,folder in [('FJ1486',parent/'cap_swivel_flexor'),('FJ1477',parent/'cap_swivel_study')]:
        screen=json.loads((folder/'screen.json').read_text());assert screen['selected_fork'] is not None and not screen['local_rotor_hardware_intersections'];C=np.array(screen['pivot_zero_display_mm']);e=np.array(screen['trunnion_axis_zero_display']);d=json.loads((muscle_folder(fid)/'design.json').read_text());z=np.array(d['basis_columns_display'])[:,2];offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));A=np.array(d['pad_centres_display_mm'][0])+offset;selected=screen['selected_fork'];delta=np.array(selected.get('bridge_delta_display_mm',np.array(selected['bridge_direction_display'])*selected.get('bridge_reach_mm',20)))
        rotor_name='flexor_cap_rotor' if fid=='FJ1486' else 'extensor_cap_rotor';key=fid+'_cap_rotor';world[key]=load(folder/(rotor_name+'.obj'))[0];rows[key]={'key':key,'file':str((folder/(rotor_name+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':fid+'_distal_pivot'}
        moving=[k for k in rows if k.startswith(fid+'_distal_') or k==fid+'_active_wire_crimp_barrel']+[key]
        for k in moving:rows[k]['owner']=fid+'_distal_pivot'
        modules[fid]={'pivot':C,'axis':e,'body_z':z,'proximal':A,'root':C+delta,'moving':moving,'folder':str(folder.relative_to(ROOT))};forks.append(load(folder/'fork.obj')[0]);voids.append(rod(C-e*30,C+e*30,1.1))
        for name in ['pin_-1','pin_1','washer_-1','washer_1']:
            key=fid+'_'+name;world[key]=load(folder/(name+'.obj'))[0];fixed_new[key]=world[key];rows[key]={'key':key,'file':str((folder/(name+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':'elbow'}
        sources += [folder/'screen.json',muscle_folder(fid)/'design.json']+[folder/(k+'.obj') for k in [rotor_name,'fork','pin_-1','pin_1','washer_-1','washer_1']]
    # Extra bony context stays in source coordinates, with the same registration
    # as the previous complete rigid arm audit.
    refs={}
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in rows:
            p=OUT/'repaired_supports'/(r['id']+'.obj');refs['reference_'+r['id']]=load(p,raw=True)[0];sources.append(p)
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    def cap_pose(fid,angle):
        m=modules[fid];R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();C=m['pivot'];e=m['axis'];z=m['body_z'];centre=E+R@(C-E);target=unit(centre-m['proximal']);target=R.T@target;target=unit(target-e*(target@e));psi=float(np.arctan2(e@np.cross(z,target),z@target));Q=R@Rotation.from_rotvec(e*psi).as_matrix();return np.c_[Q,centre-Q@C],psi
    angles=[float(q) for q in range(151)];adaptive=[angles[0]]
    def refine(lo,hi):
        change=max([hi-lo]+[abs(np.degrees(cap_pose(fid,hi)[1]-cap_pose(fid,lo)[1])) for fid in modules])
        if change>1.000001:
            mid=(lo+hi)/2;refine(lo,mid);refine(mid,hi)
        else:adaptive.append(hi)
    for lo,hi in zip(angles,angles[1:]):refine(lo,hi)
    poses=[]
    for q in adaptive:
        R=Rotation.from_rotvec(axis*np.radians(q)).as_matrix();Tf=np.c_[R,E-R@E];Ts={fid:cap_pose(fid,q)[0] for fid in modules};transformed={}
        for key,s in world.items():
            owner=rows[key]['owner'];fid=next((fid for fid in modules if owner==fid+'_distal_pivot'),None)
            transformed[key]=s.transform(Ts[fid]) if fid else s.transform(Tf) if owner=='elbow' else s
        poses.append((q,Tf,Ts,transformed))
    roots=[];ears=[]
    for anchor in j['bone_anchors']:
        q=np.array(anchor['front_spotface_mm']);roots.append(q+n*1.5);ears.append(rod(q,q+n*3,5)-rod(q-n*8,q+n*20,1.7));voids += [rod(q-n*2,q+n*30,1.7),rod(q+n*3,q+n*30,3.75)]
    all_voids=mf.Manifold.batch_boolean(voids,mf.OpType.Add);trials=[];chosen=None
    fixed_keys=[k for k in world if rows[k]['kind']!='TPU'];changed_moving=[k for m in modules.values() for k in m['moving']]
    # First check the two simultaneously moving cap packs and new fixed pins.
    # Old stable rigid pairs already passed the parent audit, but every changed
    # cap/pin pair is checked again, including adjacent pieces.
    changed=set(changed_moving)|set(fixed_new);pair_first=None
    for q,Tf,Ts,transformed in poses:
        hits=[]
        for ka,kb in itertools.combinations(fixed_keys,2):
            if ka not in changed and kb not in changed:continue
            v=intersection(transformed[ka],transformed[kb])
            if v>.001:hits.append({'pair':[ka,kb],'volume_mm3':v})
        for key in changed:
            for other,s in refs.items():
                v=intersection(transformed[key],s)
                if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
        if hits:pair_first={'elbow_deg':q,'hits':hits};break
    print('paired cap packs first limit',pair_first,'adaptive poses',len(poses),flush=True)
    if pair_first is None:
        for reach,dodge,extra in itertools.product([10,16,22,28],[-18,18,-26,26],[0,10,20]):
            exterior=[root+down*8+n*reach+side*dodge for root in roots];centre=np.mean(exterior,axis=0);pieces=forks+ears+[rod(exterior[0],exterior[1],3)];paths=[]
            for root,end in zip(roots,exterior):
                points=[root,root+down*8,root+down*8+n*reach,end];pieces.append(path(points));paths.append([p.tolist() for p in points])
            for fid,m in modules.items():
                target=m['root'];stage=centre.copy();stage[1]=target[1]+(-1 if fid=='FJ1486' else 1)*extra;stage[2]=target[2]-10;points=[centre,stage,target];pieces.append(path(points));paths.append([p.tolist() for p in points])
            frame=(mf.Manifold.batch_boolean(pieces,mf.OpType.Add)-all_voids).as_original().simplify(.001);components=len([s for s in frame.decompose() if s.volume()>.001]);first=None
            if components==1:
                for q,Tf,Ts,transformed in poses:
                    current=frame.transform(Tf);hits=[]
                    for key in fixed_keys:
                        v=intersection(current,transformed[key])
                        if v>.001:hits.append({'pair':['paired_bone_frame',key],'volume_mm3':v})
                    for key,s in refs.items():
                        v=intersection(current,s)
                        if v>.001:hits.append({'pair':['paired_bone_frame',key],'volume_mm3':v})
                    if hits:first={'elbow_deg':q,'hits':hits};break
            trial={'reach_mm':reach,'side_dodge_mm':dodge,'extra_mm':extra,'positive_components':components,'paths_display_mm':paths,'first_rigid_interference':first};trials.append(trial);print('frame trial',reach,dodge,extra,'first',first,flush=True)
            if components==1 and first is None:chosen=(frame,trial);break
    report={'scope':'real connected bone-rooted support and two rotating distal cap packs; all changed rigid pairs vs parent arm and extra bones on conditional elbow path; no TPU deformation/passive equilibrium orfull3Dshoulder audit','source_sha256':hashes,'sample_count':len(poses),'sampling':'adaptive, elbow and each prescribed passive cap rotation step <=1deg','paired_cap_pack_first_interference':pair_first,'trials':trials,'selected':chosen[1] if chosen else None,'original_bone_seats_preserved':True,'bone_centres_or_elbow_axis_changed':False,'TPU_deformation_verified':False,'installed_complete_assembly_verified':False,'frozen':False};(dest/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    if chosen:
        frame,trial=chosen;newmeta.append(write(dest/'paired_bone_frame.stl',frame));rows['paired_bone_frame']={'key':'paired_bone_frame','file':str((dest/'paired_bone_frame.obj').relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':'elbow'}
        for fid,m in modules.items():
            source=json.loads((ROOT/m['folder']/'cad_input.json').read_text());rotor_name='flexor_cap_rotor' if fid=='FJ1486' else 'extensor_cap_rotor'
            for p in source['parts']:
                if p['name']=='fork':continue
                p=dict(p);p['name']=fid+'_cap_rotor' if p['name']==rotor_name else fid+'_'+p['name'];newmeta.append(p)
        (dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':newmeta},indent=2),encoding='utf8');out=dict(a);out['instances']=list(rows.values());out.update(scope='paired passive distal cap hinges on real bone-rooted support; conditional geometry, compliant muscle unverified',frozen=False);(dest/'assembly.json').write_text(json.dumps(out,indent=2),encoding='utf8')
        state=[]
        for q,Tf,Ts,_ in poses:state.append({'elbow_deg':q,'forearm_pose':Tf.tolist(),'cap_poses':{fid:T.tolist() for fid,T in Ts.items()},'relative_cap_angles_deg':{fid:float(np.degrees(cap_pose(fid,q)[1])) for fid in modules}})
        (dest/'conditional_motion.json').write_text(json.dumps({'scope':'prescribed cap orientation, not passive equilibrium or driveresult','poses':state},indent=2),encoding='utf8');print('paired support selected; manufactured',len(rows),'newCADparts',len(newmeta),flush=True)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())

if __name__=='__main__':run()
