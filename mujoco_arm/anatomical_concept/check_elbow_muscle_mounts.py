"""M2 neutral manufactured assembly, hard reference context, and wire diagnostics."""
from build_elbow_muscle_mounts import *
from build_shaped_elbow_muscles import capsule

def run():
    data=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));solids={};boxes={};bykey={r['key']:r for r in data['instances']}
    for r in data['instances']:
        solid=load(ROOT/r['file'])[0].translate(r['translation_display_mm']);solids[r['key']]=solid;boxes[r['key']]=np.array(solid.bounding_box())
    def intersection(a,sa,ba,b,sb,bb):
        if np.any(np.maximum(ba[:3]-bb[3:],bb[:3]-ba[3:])>=0):return None
        volume=max(0,float((sa^sb).volume()))
        return {'parts':[a,b],'volume_mm3':volume} if volume>.001 else None
    hard=[];soft=[]
    for a,b in itertools.combinations(solids,2):
        hit=intersection(a,solids[a],boxes[a],b,solids[b],boxes[b])
        if hit:(soft if 'TPU' in [bykey[a]['kind'],bykey[b]['kind']] else hard).append(hit)
    context=[];support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in support['models']:
        if r['kind']!='bone_reference' or r['id']=='FJ3290':continue
        sb=load(OUT/'repaired_supports'/(r['id']+'.obj'),True)[0];bb=np.array(sb.bounding_box())
        for a,sa in solids.items():
            if bykey[a]['kind']=='TPU' or bykey[a]['owner']=='world':continue
            hit=intersection(a,sa,boxes[a],r['id'],sb,bb)
            if hit:context.append(hit)
    wires=[]
    for p in data['placements']:
        design=json.loads((MUSCLES/p['id']/'design.json').read_text(encoding='utf8'));nodes=np.array(design['wire_polyline_display_mm'])+p['translation_display_mm'];wire=mf.Manifold.batch_boolean([capsule(a,b,.3) for a,b in zip(nodes,nodes[1:])],mf.OpType.Add);wb=np.array(wire.bounding_box());hits=[]
        for key,s in solids.items():
            if bykey[key]['kind']=='TPU':continue
            hit=intersection(p['id']+'_wire',wire,wb,key,s,boxes[key])
            if hit:hits.append(hit)
        wires.append({'id':p['id'],'nominal_wire_interferences':hits,'scope':'neutral active wire polyline, not passive output or motor/conduit routing'})
    report={'scope':'all100 manufactured neutral instances checked; new compliant operation and drive NOT verified','manufactured_instance_count':len(solids),'all_manufactured_rigid_pair_interferences':hard,'TPU_contact_review_flags':soft,'TPU_policy':'soft intersections require deformation/friction review; neither automatic pass nor hard-failure classification','rigid_vs_fixed_reference_interferences':context,'nominal_active_wires':wires,'volume_threshold_mm3':.001,'neutral_rigid_geometry_clear':not hard and not context,'frozen':False}
    report['scope']=f'all{len(solids)} manufactured neutral instances; explicit hard/soft classification, no compliant operation or drive proof'
    (DEST/'neutral_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('Neutral',len(solids),'hard pairs',len(hard),'context',len(context),'soft flags',len(soft),'wire hits',[len(r['nominal_wire_interferences']) for r in wires],flush=True)
    for hit in hard[:12]:print(hit,flush=True)

if __name__=='__main__':run()
