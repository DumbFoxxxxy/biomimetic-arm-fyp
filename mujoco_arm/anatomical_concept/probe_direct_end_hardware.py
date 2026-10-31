"""Rigid-space screen before constructing direct forearm/TPU-end adapters.

Conditionally carries distal caps with forearm; no adapter exists yet. Does not
prescribe, deform, or approve the TPU lattice. Includes all hard adjacent bones.
"""
from build_forearm_receivers import ROOT,OUT,BASE,load,np,mf,intersection,Rotation
import json,hashlib
DEST=OUT/'direct_muscle_interfaces'
M2=OUT/'elbow_muscles/mounted_elbow_candidate'

def run():
    a=json.loads((M2/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);distals={k:s for k,s in world.items() if any(k.startswith(fid+'_distal_') or k==fid+'_active_wire_crimp_barrel' for fid in ['FJ1486','FJ1477'])};others={k:s for k,s in world.items() if k not in distals and rows[k]['kind']!='TPU'}
    supports=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in supports['models']:
        if r['kind']=='bone_reference' and r['id'] not in others:
            k='reference_'+r['id'];others[k]=load(OUT/'repaired_supports'/(r['id']+'.obj'),raw=True)[0];rows[k]={'owner':'world'}
    reports=[]
    for fid in ['FJ1486','FJ1477']:
        first=None;poses=[];keys=[k for k in distals if k.startswith(fid+'_')]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];hits=[]
            for key in keys:
                moved=distals[key].transform(T)
                for k,s in others.items():
                    body=s.transform(T) if rows[k]['owner']=='elbow' else s;v=intersection(moved,body)
                    if v>.001:hits.append({'pair':[key,k],'volume_mm3':v})
            if hits and first is None:first={'elbow_angle_deg':angle,'hits':hits}
            poses.append({'angle_deg':angle,'hits':hits})
        reports.append({'muscle':fid,'rigid_pack_components':len(keys),'first_rigid_collision':first,'poses':poses})
        print(fid,'direct end hard-packing first limit',first,flush=True)
    (DEST/'direct_end_hardware_probe.json').write_text(json.dumps({'scope':'conditional distal cap/hardware fixed to forearm, adapter not constructed. TPU omitted from hard test because deformation is unknown; no muscle movement/drive pass.','modules':reports,'step_deg':1,'source_assembly_sha256':hashlib.sha256((M2/'assembly.json').read_bytes()).hexdigest(),'drive_verified':False,'frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
