"""Neutral complete geometry and independent return of changed rigid paths.

Unchanged rigid pairs reuse the parent's source-guarded elbow path evidence.
No flexible equilibrium, actual passive joint angle or drive range claim.
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

def run():
    dest=DEST/'paired_swivel_support'
    if '--keepers' in sys.argv:dest=dest/'captive_keepers'
    motion_path=(dest.parent if '--keepers' in sys.argv else dest)/'conditional_motion.json'
    a=json.loads((dest/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};inputs=[dest/'assembly.json',dest/'design.json',motion_path]+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};hard=[];soft=[]
    for (ka,sa),(kb,sb) in itertools.combinations(world.items(),2):
        v=intersection(sa,sb)
        if v>.001:(soft if 'TPU' in [rows[ka]['kind'],rows[kb]['kind']] else hard).append({'pair':[ka,kb],'volume_mm3':v})
    pose_data=json.loads(motion_path.read_text())['poses'];q=np.array([p['elbow_deg'] for p in pose_data]);steps=[float(abs(np.diff(q)).max())]+[float(abs(np.diff([p['relative_cap_angles_deg'][fid] for p in pose_data])).max()) for fid in ['FJ1486','FJ1477']];assert max(steps)<=1.000001
    rigid=[k for k,r in rows.items() if r['kind']!='TPU' and not r.get('nominal_packing_only')];changed={k for k in rigid if k=='paired_bone_frame' or k.endswith('_cap_rotor') or any(k.startswith(fid+'_'+prefix) for fid in ['FJ1486','FJ1477'] for prefix in ['distal_','pin_','washer_','keeper'])};pairs=[(ka,kb) for ka,kb in itertools.combinations(rigid,2) if ka in changed or kb in changed];refs={}
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in rows:
            p=OUT/'repaired_supports'/(r['id']+'.obj');refs['reference_'+r['id']]=load(p,raw=True)[0];hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    parent=json.loads((DEST/'bone_motion_audit.json').read_text());assert parent['complete_prescribed_suite'];assert next(x for x in parent['actions'] if x['probe']=='elbow_flexion')['continuous_prefix_fraction']==1
    for name,h in parent['input_evidence_sha256'].items():assert hashlib.sha256((DEST/name).read_bytes()).hexdigest()==h
    first=None
    for pose in reversed(pose_data):
        transformed={}
        for key in rigid:
            owner=rows[key]['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else np.array(pose['forearm_pose']) if owner=='elbow' else None;transformed[key]=world[key].transform(T) if T is not None else world[key]
        hits=[]
        for ka,kb in pairs:
            v=intersection(transformed[ka],transformed[kb])
            if v>.001:hits.append({'pair':[ka,kb],'volume_mm3':v})
        for key in changed:
            for other,s in refs.items():
                v=intersection(transformed[key],s)
                if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
        if hits:first={'elbow_deg':pose['elbow_deg'],'hits':hits};break
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    result={'scope':f'all{len(rows)} nominal neutral manufactured instances, independently recomputed return of all changed rigid pairs and extra bones; unchanged parent rigid-pair sourceproof reused; no passive equilibrium, flexible motion oractuation','neutral_rigid_interferences':hard,'neutral_TPU_contact_flags':soft,'motion_attached_rigid_count':len(rigid),'changed_rigid_count':len(changed),'changed_pair_count':len(pairs),'return_samples':len(pose_data),'maximum_elbow_and_cap_step_deg':max(steps),'first_return_interference':first,'parent_rigid_source_evidence_sha256':hashlib.sha256((DEST/'bone_motion_audit.json').read_bytes()).hexdigest(),'source_sha256':hashes,'complete_TPU_assembly_or_drive_verified':False,'frozen':False};(dest/'geometry_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('paired neutral hard',len(hard),'TPU flags',len(soft),'return first',first,'samples',len(pose_data),flush=True)

if __name__=='__main__':run()
