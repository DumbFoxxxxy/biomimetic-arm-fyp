"""Staged keeper/forearm subassembly paths, no complete compliant arm claim."""

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
    dest=DEST/'paired_swivel_support/captive_keepers';a=json.loads((dest/'assembly.json').read_text());design=json.loads((dest/'design.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};inputs=[dest/'assembly.json',dest/'design.json']+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};results=[]
    def probe(label,moving,fixed,direction,length):
        first=None;soft=[]
        for distance in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
            hits=[]
            for key in moving:
                s=world[key].translate(direction*distance)
                for other in fixed:
                    v=intersection(s,world[other])
                    if v>.001:
                        h={'pair':[key,other],'volume_mm3':v}
                        if 'TPU' in [rows[key]['kind'],rows[other]['kind']]:soft.append({'distance_mm':float(distance),**h})
                        else:hits.append(h)
            if hits:first={'distance_mm':float(distance),'hits':hits};break
        result={'label':label,'moving':moving,'fixed':fixed,'direction_display':direction.tolist(),'approach_mm':length,'rigid_clear':first is None,'first_rigid_interference':first,'soft_flags':soft,'final5mm_step_mm':.1,'other_step_mm':1};results.append(result);print(label,'clear',first is None,'first',first,flush=True);return result
    modules=['FJ1486','FJ1477'];cap_pack=[k for k in world if any(k.startswith(fid+'_distal_') or k==fid+'_cap_rotor' for fid in modules)];pins=[fid+'_pin_'+str(s) for fid in modules for s in [-1,1]];keepers=[fid+'_keeper' for fid in modules];clamps=[k for k in world if any(k.startswith(fid+'_keeper_') for fid in modules)];base=['paired_bone_frame']+cap_pack
    for fid,folder in [('FJ1486','cap_swivel_flexor'),('FJ1477','cap_swivel_study')]:
        screen=json.loads((DEST/folder/'screen.json').read_text());e=np.array(screen['trunnion_axis_zero_display'])
        for s in [-1,1]:
            key=fid+'_pin_'+str(s);probe(key+'_bench_entry',[key],base+[k for k in pins if k!=key],e*s,12)
    for feature in design['features']:
        fid=feature['muscle'];direction=np.array(feature['keeper_bridge_direction_display']);axis=np.array(feature['keeper_bolt_axis_display']);key=fid+'_keeper'
        back=[fid+'_keeper_nut',fid+'_keeper_back_washer'];front=[fid+'_keeper_bolt',fid+'_keeper_front_washer'];fixed=base+pins+[k for k in clamps if k not in back+front];screen=json.loads((DEST/('cap_swivel_flexor' if fid=='FJ1486' else 'cap_swivel_study')/'screen.json').read_text());e=np.array(screen['trunnion_axis_zero_display']);normal=unit(np.cross(axis,e));trials=[]
        for vector,name in [(e,'axis_plus'),(-e,'axis_minus'),(normal,'normal_plus'),(-normal,'normal_minus')]:
            trial=probe(fid+'_rear_nut_before_keeper_'+name,back,fixed,vector,20);trials.append(trial)
            if trial['rigid_clear']:break
        rear_all=[k for k in clamps if k.endswith(('_keeper_nut','_keeper_back_washer'))]
        probe(key+'_side_entry_after_rear_nuts',[key],base+pins+rear_all+[k for k in keepers if k!=key],direction,25)
        probe(fid+'_front_clamp_after_rear_nut',front,base+pins+keepers+[k for k in clamps if k not in front],axis,25)
    unit_keys=base+pins+keepers+clamps;j=json.loads((BASE/'design.json').read_text());n=np.array(j['bone_anchors'][0]['bolt_axis']);down=np.array(j['zero_forearm_direction']);mount_trials=[]
    # Bare forearm subassembly: humerus, motor, TPU and elbow joint pin are not
    # already attached. This order cannot certify the entire arm installation.
    for direction,length,label in [(n,60,'normal_front'),(-n,60,'normal_back'),(down,320,'from_distal_forearm_end'),(-down,120,'from_proximal_forearm_end')]:
        trial=probe('rigid_unit_to_bare_forearm_'+label,unit_keys,['FJ3391'],direction,length);mount_trials.append(trial)
        if trial['rigid_clear']:break
    if any(r['rigid_clear'] for r in mount_trials):
        fixed=unit_keys+['FJ3391']
        for i in [0,1]:
            keys=[f'bone_{i}_bolt',f'bone_{i}_front_washer'];probe(f'bone_{i}_front_after_unit',keys,fixed,n,45);fixed+=keys;keys=[f'bone_{i}_nut',f'bone_{i}_back_washer'];probe(f'bone_{i}_rear_after_unit',keys,fixed,-n,35);fixed+=keys
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    report={'scope':'local rigid cap/pin/keeper bench fitting, then rigid unit to bare coupledforearm and bonebolts; humerus/TPU/steel/actuator/elbow-joint installation, tools, thread rotation andpreload NOT verified','paths':results,'bare_forearm_mount_found':any(r['rigid_clear'] for r in mount_trials),'source_sha256':hashes,'complete_arm_assembly_verified':False,'frozen':False};(dest/'assembly_paths.json').write_text(json.dumps(report,indent=2),encoding='utf8')

if __name__=='__main__':run()
