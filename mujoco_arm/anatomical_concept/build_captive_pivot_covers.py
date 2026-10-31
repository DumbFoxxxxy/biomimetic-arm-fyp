"""Actual custom helical retention covers around the four split shoulder pins.

Not an ISO-thread designation. Geometry establishes a captive pin head and
remaining blind-bearing engagement; strength, printing and loosening untested.
Original M10 files remain unchanged in their parent directory.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from project_runtime import WorkspacePath as Path
from scipy.spatial.transform import Rotation

def helical_ridge(profile,start=18.65,end=20.75,pitch=1.):
    theta=np.linspace(0,(end-start)*2*np.pi/pitch,int(np.ceil((end-start)/pitch*96))+1);vertices=[]
    for t in theta:
        for r,dz in profile:vertices.append([r*np.cos(t),r*np.sin(t),start+t*pitch/(2*np.pi)+dz])
    vertices=np.array(vertices);faces=[];n=len(profile)
    for i in range(len(theta)-1):
        for j in range(n):a=i*n+j;b=i*n+(j+1)%n;c=(i+1)*n+(j+1)%n;d=(i+1)*n+j;faces.extend([[a,b,c],[a,c,d]])
    for j in range(1,n-1):faces.append([0,j+1,j]);faces.append([(len(theta)-1)*n,(len(theta)-1)*n+j,(len(theta)-1)*n+j+1])
    faces=np.array(faces,dtype=np.uint64);s=mf.Manifold(mf.Mesh64(vertices.astype(np.float64),faces))
    if s.volume()<0:s=mf.Manifold(mf.Mesh64(vertices.astype(np.float64),faces[:,::-1].copy()))
    assert s.status()==mf.Error.NoError and s.volume()>0
    return s

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';dest=parent/'captive_covers';dest.mkdir(exist_ok=True);a=json.loads((parent/'assembly.json').read_text());rows={r['key']:dict(r) for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};frame=world['paired_bone_frame'];baseframe=frame;screens={fid:json.loads((DEST/folder/'screen.json').read_text()) for fid,folder in [('FJ1486','cap_swivel_flexor'),('FJ1477','cap_swivel_study')]};cover_local=(rod([0,0,18.35],[0,0,22.4],5,64)-rod([0,0,18],[0,0,21.35],3.6,64)-rod([0,0,21.25],[0,0,22.5],1.2,48));male_local=(rod([0,0,18.1],[0,0,21.2],3.5,64)-rod([0,0,18],[0,0,21.3],2.7,64))+helical_ridge([(3.45,-.32),(4.,-.125),(4.,.125),(3.45,.32)])
    groove=helical_ridge([(3.4,-.42),(4.12,-.225),(4.12,.225),(3.4,.42)]);cover_local=(cover_local-groove).as_original().simplify(.0001)
    for s in [-1,1]:cover_local=cover_local-rod([s*3.2,0,21.2],[s*3.2,0,22.5],.6,32)
    male_local=male_local.as_original().simplify(.0001);covers={};features=[];sources=[parent/'assembly.json',parent/'conditional_motion.json']+[ROOT/r['file'] for r in rows.values()]
    for fid,screen in screens.items():
        C=np.array(screen['pivot_zero_display_mm']);e=np.array(screen['trunnion_axis_zero_display']);d=json.loads((muscle_folder(fid)/'design.json').read_text());z=np.array(d['basis_columns_display'])[:,2]
        for sign in [-1,1]:
            v=e*sign;F=np.column_stack([z,np.cross(v,z),v]);T=np.c_[F,C];male=male_local.transform(T);cap=cover_local.transform(T);frame=frame+male;key=fid+'_retention_cover_'+str(sign);covers[key]=cap
            features.append({'key':key,'muscle':fid,'sign':sign,'centre_display_mm':C.tolist(),'axis_display':v.tolist(),'basis_columns':F.tolist(),'male_bounds_outward_mm':[18.1,21.2],'cap_bounds_outward_mm':[18.35,22.4],'pitch_mm':1.,'nominal_effective_thread_turns':2.1,'head_outward_float_mm':1.35,'minimum_pin_tip_to_stub_end_engagement_mm':1.95})
        sources += [DEST/('cap_swivel_flexor' if fid=='FJ1486' else 'cap_swivel_study')/'screen.json',muscle_folder(fid)/'design.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};frame=frame.as_original().simplify(.0001);components=len([s for s in frame.decompose() if s.volume()>.001]);assert components==1;parts={'paired_bone_frame':frame,**covers};world.update(parts);hard=[];soft=[]
    for key,s in parts.items():
        for other,body in world.items():
            if other==key or other in parts and other<key:continue
            volume=intersection(s,body)
            if volume>.001:(soft if rows.get(other,{}).get('kind')=='TPU' else hard).append({'pair':[key,other],'volume_mm3':volume})
    captures=[]
    for feature in features:
        fid=feature['muscle'];sign=feature['sign'];v=np.array(feature['axis_display']);pin=world[fid+'_pin_'+str(sign)];cover=covers[feature['key']];washer=world[fid+'_washer_'+str(sign)];captures.append({'key':feature['key'],'pin_rest_cover_overlap_mm3':intersection(pin,cover),'pin_outward_1_34mm_cover_overlap_mm3':intersection(pin.translate(v*1.34),cover),'pin_outward_1_36mm_cover_overlap_mm3':intersection(pin.translate(v*1.36),cover),'pin_inward_0_01mm_washer_overlap_mm3':intersection(pin.translate(-v*.01),washer),'remaining_blind_bearing_engagement_at_stop_mm':1.95})
    print('captive covers neutral hard',len(hard),'TPU',len(soft),'frame components',components,flush=True)
    motion=json.loads((parent/'conditional_motion.json').read_text())['poses'];first=None
    # All new/revised fixed pieces vs all existing hard pieces, including both
    # rotating cap packs. Recheck the extra bone references as in the parent.
    refs={}
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in rows:
            p=OUT/'repaired_supports'/(r['id']+'.obj');refs['reference_'+r['id']]=load(p,raw=True)[0];hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for pose in motion:
        Tf=np.array(pose['forearm_pose']);hits=[]
        for key,s in parts.items():
            current=s.transform(Tf)
            for other,body in world.items():
                if other==key or other in parts or rows.get(other,{}).get('kind')=='TPU':continue
                owner=rows[other]['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;other_s=body.transform(T) if T is not None else body;volume=intersection(current,other_s)
                if volume>.001:hits.append({'pair':[key,other],'volume_mm3':volume})
            for other,body in refs.items():
                volume=intersection(current,body)
                if volume>.001:hits.append({'pair':[key,other],'volume_mm3':volume})
        if hits:first={'elbow_deg':pose['elbow_deg'],'hits':hits};break
    report={'scope':'custom actual helical collars and covers around existing four pins; nominal capture geometry plus changed rigid swept clearance, not strength/loosening/print qualification','source_sha256':hashes,'thread_profile':'custom8mmcrest diameter/1mmpitch, notISOdesignation; male3.5..4.0mmR, female bore3.6 andgroove4.12mmR','features':features,'frame_components':components,'neutral_rigid_interferences':hard,'neutral_TPU_review_flags':soft,'capture_witnesses':captures,'conditional_samples':len(motion),'first_conditional_interference':first,'thread_screw_insertion_verified':False,'vibration_loosening_verified':False,'physical_retention_verified':False,'frozen':False};(dest/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('captive swept first',first,'capture stops',[r['pin_outward_1_36mm_cover_overlap_mm3'] for r in captures],flush=True)
    if not hard and not soft and first is None:
        metas=[]
        for key,s in parts.items():
            metas.append(write(dest/(key+'.stl'),s));rows[key]={'key':key,'file':str((dest/(key+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':'elbow'}
        out=dict(a);out['instances']=list(rows.values());out.update(scope='M10 retention-cover candidate, actual customhelical capture, no printing/strength/loosening proof',frozen=False);(dest/'assembly.json').write_text(json.dumps(out,indent=2),encoding='utf8');(dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':metas},indent=2),encoding='utf8');print('captive candidate',len(rows),'instances,5revised/newCADparts',flush=True)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())

if __name__=='__main__':run()
