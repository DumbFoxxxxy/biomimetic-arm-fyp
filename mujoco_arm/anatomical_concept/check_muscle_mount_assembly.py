"""Finite mounting sequence for M2, distinct from muscle operation/deformation."""
from build_elbow_muscle_mounts import *

def run():
    data=json.loads((DEST/'assembly.json').read_text());solids={};rows={r['key']:r for r in data['instances']}
    for r in data['instances']:solids[r['key']]=load(ROOT/r['file'])[0].translate(r['translation_display_mm'])
    def probe(label,moving,fixed,direction,distance):
        samples=[]
        for offset in range(distance,-1,-1):
            hits=[]
            for a in moving:
                moved=solids[a].translate(direction*offset);box=np.array(moved.bounding_box())
                for b in fixed:
                    other=solids[b];bb=np.array(other.bounding_box())
                    if np.any(np.maximum(box[:3]-bb[3:],bb[:3]-box[3:])>=0):continue
                    volume=max(0,float((moved^other).volume()))
                    if volume>.001:hits.append({'parts':[a,b],'volume_mm3':volume,'classification':'soft_review' if 'TPU' in [rows[a]['kind'],rows[b]['kind']] else 'rigid'})
            samples.append({'offset_mm':offset,'hits':hits})
        return {'label':label,'moving':moving,'direction':direction.tolist(),'step_mm':1,'samples':samples,'rigid_path_clear':not any(h['classification']=='rigid' for s in samples for h in s['hits']),'soft_contact_flags':sum(h['classification']=='soft_review' for s in samples for h in s['hits'])}
    results=[];skeleton=[k for k in solids if not k.startswith(('FJ1486_','FJ1477_'))];installed=[]
    for p in data['placements']:
        fid=p['id'];n=np.array(p['bone_anchors'][0]['bolt_axis']);design=json.loads((MUSCLES/fid/'design.json').read_text());z=np.array(design['basis_columns_display'])[:,2]
        adapter=fid+'_bone_adapter'
        cap_fasteners=[fid+'_proximal_'+suffix for suffix in ['bolt_left','bolt_right','head_washer_left','head_washer_right']]
        module=[k for k in solids if k.startswith(fid+'_') and '_bone_' not in k and k not in cap_fasteners]
        # Assemble the module/adapter off the bone before transverse bone bolts.
        results.append(probe(fid+'_module_side_entry_to_adapter_on_bench',module,[adapter],n,25))
        for side in ['left','right']:
            moving=[fid+'_proximal_bolt_'+side,fid+'_proximal_head_washer_'+side];results.append(probe(fid+'_shared_cap_bolt_on_bench_'+side,moving,[k for k in module+cap_fasteners+[adapter] if k not in moving],-z,25))
        unit=module+cap_fasteners+[adapter]
        results.append(probe(fid+'_preassembled_unit_to_bone',unit,skeleton+installed,n,30));installed+=unit
    # Bone bolts enter after both units have seated. This does not presume that
    # uncharacterized TPU intersections are automatically acceptable.
    for p in data['placements']:
        fid=p['id'];n=np.array(p['bone_anchors'][0]['bolt_axis'])
        for i,a in enumerate(p['bone_anchors']):
            prefix=fid+'_bone_'+str(i);moving=[prefix+'_bolt',prefix+'_front_washer'];fixed=[k for k in solids if k not in moving and k not in [prefix+'_nut',prefix+'_back_washer']]
            results.append(probe(prefix+'_bolt_front_after_units',moving,fixed,n,30));moving=[prefix+'_nut',prefix+'_back_washer'];results.append(probe(prefix+'_nut_back_after_units',moving,[k for k in solids if k not in moving],-n,30))
    report={'scope':'sampled nominal assembly, unloaded module form manually held during fitting; not compliant operation, threads/preload, tool access or forces','sequence':['assemble module and adapter on bench, transverse bone bolts absent','insert shared M3x16 cap fasteners on bench','seat complete module/adapter units onto bone','insert transverse bone bolts and rear washer/nuts'],'paths':results,'all_rigid_paths_clear':all(r['rigid_path_clear'] for r in results),'soft_flags':sum(r['soft_contact_flags'] for r in results),'complete_assembly_verified':False,'frozen':False}
    (DEST/'mount_assembly_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('mount fitting',sum(r['rigid_path_clear'] for r in results),'/',len(results),'soft flags',report['soft_flags']);print('failed',[r['label'] for r in results if not r['rigid_path_clear']])

if __name__=='__main__':run()
