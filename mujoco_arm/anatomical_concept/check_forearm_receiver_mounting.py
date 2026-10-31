"""Changed lever and receiving-unit mounting paths in full neutral context."""
from build_forearm_receivers import *

MOUNT=DEST/'rephased_anchors'

def run():
    a=json.loads((PARENT/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};source_paths=[PARENT/'assembly.json',BASE/'design.json',MOUNT/'mounted_direction_screen.json',MOUNT/'lever_receiving_spotfaces.obj']+[ROOT/r['file'] for r in rows.values()]+list(MOUNT.glob('*_neutral.obj'));source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if k not in ['side_lever','flexor_tendon_stop','extensor_tendon_stop']};world['new_lever']=load(MOUNT/'lever_receiving_spotfaces.obj')[0];probe_data=json.loads((MOUNT/'mounted_direction_screen.json').read_text());definition=json.loads((BASE/'design.json').read_text());axis=np.array(definition['joint_axis']);n=-axis;down=np.array(definition['zero_forearm_direction']);F=np.stack([down,np.cross(n,down),n],axis=1)
    for p in MOUNT.glob('*_neutral.obj'):world[p.stem.removesuffix('_neutral')]=load(p)[0]
    results=[]
    def probe(label,moving,fixed,direction,length):
        first=None;soft=[]
        for distance in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
            h=[]
            for key in moving:
                s=world[key].translate(direction*distance)
                for other in fixed:
                    v=intersection(s,world[other])
                    if v>.001:
                        hit={'pair':[key,other],'volume_mm3':v}
                        if other in rows and rows[other]['kind']=='TPU':soft.append({'distance_mm':float(distance),**hit})
                        else:h.append(hit)
            if h and first is None:first={'distance_mm':float(distance),'hits':h}
        r={'label':label,'moving':moving,'direction':direction.tolist(),'first_limit':first,'soft_flags':soft,'rigid_clear':first is None,'final5mm_step_mm':.1,'other_step_mm':1};results.append(r);print(label,'clear',r['rigid_clear'],'soft',len(soft),flush=True);return r
    fasteners=[f'bone_{i}_{suffix}' for i in [0,1] for suffix in ['bolt','front_washer','nut','back_washer']];newnames=[k for k in world if k not in rows];fixed=[k for k in world if k not in newnames+fasteners]
    probe('rephased_lever_before_bolts_and_receivers',['new_lever'],fixed,n,60);fixed.append('new_lever')
    for i in [0,1]:
        prefix='bone_'+str(i);pair=[prefix+'_bolt',prefix+'_front_washer'];probe(prefix+'_front',pair,fixed,n,45);fixed+=pair
        pair=[prefix+'_nut',prefix+'_back_washer'];probe(prefix+'_rear',pair,fixed,-n,35);fixed+=pair
    for fid in ['FJ1486','FJ1477']:
        pair=[fid+'_yaw_base_washer',fid+'_yaw_sleeve'];probe(fid+'_sleeve_on_lever',pair,fixed,n,25);fixed+=pair
        group=[fid+'_'+k for k in YAW_PARTS];probe(fid+'_assembled_yoke_onto_yaw_sleeve',group,fixed,n,30);fixed+=group
        pair=[fid+'_yaw_back_washer',fid+'_yaw_nut'];options=[]
        for label,direction in [('radial',unit(np.array(probe_data['candidate_eyelets_display_mm']['flexor' if fid=='FJ1486' else 'extensor'])-np.array(definition['joint_centre_display_mm'])-n*56)),('down',down),('up',-down),('side',F[:,1]),('other_side',-F[:,1]),('axial',-n)]:
            r=probe(fid+'_loose_rear_nut_'+label,pair,fixed,direction,20);results.pop();options.append(r)
            if r['rigid_clear'] and not r['soft_flags']:break
        chosen=min(options,key=lambda r:(not r['rigid_clear'],len(r['soft_flags'])));chosen['tested_approaches']=[{k:r[k] for k in ['label','rigid_clear','first_limit']} for r in options];results.append(chosen);fixed+=pair
        pair=[fid+'_yaw_M3',fid+'_yaw_front_washer'];probe(fid+'_yaw_bolt_after_loose_nut',pair,fixed,n,30);fixed+=pair
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in source_hashes.items())
    result={'scope':'13 changed rigid mounting paths, neutral source muscle forms held; not complete compliant assembly, tool access, thread/preload or force proof','paths':results,'all_clear':all(r['rigid_clear'] for r in results),'soft_flag_count':sum(len(r['soft_flags']) for r in results),'source_sha256':source_hashes,'frozen':False};(MOUNT/'mounting_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8')

if __name__=='__main__':run()
