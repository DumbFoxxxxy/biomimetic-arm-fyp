"""Recheck bare elbow assembly with the real fused rotor/pocket replacement."""
from check_rotor_pocket import ROOT,DEST
from build_elbow_muscle_mounts import load,OUT,np,unit
import json,hashlib

def run():
    base=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
    data=json.loads((base/'assembly.json').read_text(encoding='utf8'))
    previous=json.loads((base/'assembly_and_retention_checks.json').read_text(encoding='utf8'))
    parts={n:load(base/(n+'.obj'))[0] for n in data['components']}
    parts['humerus_GH_elbow_fork']=load(DEST/'humerus_rotor_pocket.obj')[0]
    parts['forearm_elbow_hub']=load(DEST/'forearm_with_rotor.obj')[0]
    axis=np.array(data['axis_source']);direction=np.array(previous['chosen_cartridge_path']['direction'])
    assert previous['chosen_cartridge_path']['elbow_flexion_deg']==0
    def probe(label,moving,fixed,vector,offsets):
        samples=[]
        for distance in offsets:
            hits=[]
            for a in moving:
                moved=parts[a].translate(vector*distance)
                for b in fixed:
                    volume=max(0,float((moved^parts[b]).volume()))
                    if volume>.001:hits.append({'pair':[a,b],'volume_mm3':volume})
            samples.append({'offset_mm':float(distance),'hits':hits})
        return {'label':label,'moving':moving,'direction':vector.tolist(),'samples':samples,'clear':not any(s['hits'] for s in samples),'first_limit':next((s for s in samples if s['hits']),None)}
    hum='humerus_GH_elbow_fork';fore=['forearm_elbow_hub','elbow_bearing_bush']
    results=[probe('forearm_with_fused_rotor_and_bush_before_sleeves',fore,[hum],direction,np.r_[np.arange(45,5,-1),np.linspace(5,0,51)])]
    installed=[hum]+fore
    for name,sign in [('elbow_spacer_left',-1),('elbow_spacer_right',1)]:
        results.append(probe(name+'_axial_install',[name],installed,axis*sign,range(85,-1,-1)));installed.append(name)
    results.append(probe('pin_and_head_washer',['elbow_shoulder_pin','elbow_washer_left'],installed,-axis,range(int(np.ceil(data['custom_pin_shank_length_mm']+10)),-1,-1)))
    result={'scope':'bare nine-part elbow sequence before muscle modules installed; no threads/preload/tool-access/load proof','paths':results,'last_5mm_forearm_step_mm':.1,'other_step_mm':1,'sequence_clear':all(r['clear'] for r in results),'geometry_sha256':{n:hashlib.sha256((DEST/n).read_bytes()).hexdigest() for n in ['humerus_rotor_pocket.obj','forearm_with_rotor.obj']},'full_arm_assembly_verified':False,'frozen':False}
    (DEST/'joint_assembly_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('bare rotor elbow assembly',result['sequence_clear'],[(r['label'],r['first_limit']) for r in results if not r['clear']],flush=True)
if __name__=='__main__':run()
