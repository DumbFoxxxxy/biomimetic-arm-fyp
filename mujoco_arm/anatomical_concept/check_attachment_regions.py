"""Check numerical distance/search logic and atlas-to-model coverage."""
from extract_attachment_regions import *

def run():
    triangle=np.array([[[0.,0,0],[1.,0,0],[0.,1,0]]]);examples=[([.2,.2,1],1),([2,0,0],1),([-1,-1,0],np.sqrt(2)),([.75,.75,0],np.sqrt(.125)),([0,0,0],0)]
    for p,expected in examples:assert abs(closest(np.array(p),triangle)[1]-expected)<1e-9
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));report=json.loads((OUT/'attachment_regions/atlas_and_regions.json').read_text(encoding='utf8'));assert {r['id'] for r in report['muscles']}=={r['id'] for r in inv['models'] if r['kind']=='muscle'}
    src=next(r for r in inv['models'] if r['id']=='FJ3384');v,f=read_obj(str(ROOT/src['source']));v,f=weld_stl(v[f]);rng=np.random.default_rng(17);points=v[rng.choice(len(v),24,replace=False)]+rng.normal(0,3,(24,3));project=projector(v,f);_,fast,_=project(points);slow=np.array([closest(p,v[f])[1] for p in points]);error=float(np.max(abs(fast-slow)));assert error<1e-8
    assert all(not r['frozen'] and not r['fiber_direction_verified'] for r in report['muscles'])
    sha=all(hashlib.sha256((ROOT/r['source']).read_bytes()).hexdigest()==r['sha256'] for r in inv['models']);assert sha
    for row in report['muscles']:
        for region in row['regions']:
            if 'data' in region:
                data=np.load(OUT/'attachment_regions'/region['data']);assert len(data['distances_mm'])==region['selected_face_count_2mm'];assert np.all(data['distances_mm']<=2+1e-8)
    result={'muscle_count':len(report['muscles']),'coverage_complete_for_supplied_models':True,'distance_known_cases':len(examples),'bounded_search_vs_full_triangle_search_samples':24,'maximum_distance_discrepancy_mm':error,'all_source_hashes_unchanged':sha,'all_regions_unfrozen':True,'scope':'distance extraction integrity only; not anatomical landmark confirmation, engineering assembly or drive validation'};(OUT/'attachment_regions/checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(result)

if __name__=='__main__':run()
