"""Finite real installation paths of M4 lever/bolts; no soft assembly claim."""
from build_elbow_side_lever import *
import sys
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15';folder=DEST/variant
    design=json.loads((folder/'design.json').read_text(encoding='utf8'));assembly=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'));parts={};soft=set()
    for r in assembly['instances']:
        if r['key']=='FJ3391':continue
        parts[r['key']]=load(ROOT/r['file'])[0].translate(r['translation_display_mm'])
        if r['kind']=='TPU':soft.add(r['key'])
    new={key:load(folder/(key+'.obj'))[0] for key in design['parts']};parts.update(new)
    def probe(label,moving,fixed,direction,offsets):
        samples=[]
        for distance in offsets:
            hits=[]
            for a in moving:
                moved=parts[a].translate(direction*distance);aa=np.array(moved.bounding_box())
                for b in fixed:
                    bb=np.array(parts[b].bounding_box())
                    if np.any(np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])>=0):continue
                    volume=max(0,float((moved^parts[b]).volume()))
                    if volume>.001:hits.append({'pair':[a,b],'volume_mm3':volume,'kind':'soft_review' if a in soft or b in soft else 'rigid'})
            samples.append({'offset_mm':float(distance),'hits':hits})
        return {'label':label,'moving':moving,'direction':direction.tolist(),'samples':samples,'rigid_clear':not any(h['kind']=='rigid' for r in samples for h in r['hits']),'soft_flags':sum(h['kind']=='soft_review' for r in samples for h in r['hits']),'first_rigid_limit':next((r for r in samples if any(h['kind']=='rigid' for h in r['hits'])),None)}
    skeleton=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection';old=json.loads((skeleton/'assembly_and_retention_checks.json').read_text(encoding='utf8'));direction=np.array(old['chosen_cartridge_path']['direction']);assert old['chosen_cartridge_path']['elbow_flexion_deg']==0
    results=[probe('spotfaced_forearm_plus_bush_before_pin_sleeves_and_modules',['forearm_spotfaced','elbow_bearing_bush'],['FJ3368'],direction,np.r_[np.arange(45,5,-1),np.linspace(5,0,51)])]
    n=np.array(design['bone_anchors'][0]['bolt_axis']);new_names=set(new);base=[k for k in parts if k not in new_names]+['forearm_spotfaced']
    results.append(probe('removable_lever_radial_entry_before_bolts',['side_lever'],base,n,np.r_[np.arange(60,5,-1),np.linspace(5,0,51)]))
    fitted=base+['side_lever']
    for i in [0,1]:
        name='bone_'+str(i);pair=[name+'_bolt',name+'_front_washer']
        results.append(probe(name+'_bolt_and_front_washer',pair,fitted,n,range(45,-1,-1)));fitted+=pair
        pair=[name+'_nut',name+'_back_washer']
        results.append(probe(name+'_back_washer_and_nut',pair,fitted,-n,range(35,-1,-1)));fitted+=pair
    far=np.array(design['joint_axis'])*design['axis_side']
    for label in ['flexor','extensor']:
        key=label+'_tendon_stop';results.append(probe(key+'_rear_entry',[key],fitted,far,range(15,-1,-1)));fitted.append(key)
    report={'variant':variant,'scope':'finite geometric paths: bare forearm first, then detachable lever and nominal fasteners. Muscle modules held in existing neutral arrangement for later stages. No threads/preload/TPU handling/tool-access proof.','paths':results,'all_rigid_paths_clear':all(r['rigid_clear'] for r in results),'soft_flags':sum(r['soft_flags'] for r in results),'final5mm_step_mm':.1,'other_step_mm':1,'source_mesh_sha256':{key:hashlib.sha256((folder/(key+'.obj')).read_bytes()).hexdigest() for key in new},'complete_arm_assembly_verified':False,'frozen':False}
    (folder/'assembly_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('M4',variant,'rigid clear',sum(r['rigid_clear'] for r in results),'/',len(results),'soft',report['soft_flags'],flush=True)
    for r in results:
        if not r['rigid_clear']:print(r['label'],r['first_rigid_limit'],flush=True)
if __name__=='__main__':run()
