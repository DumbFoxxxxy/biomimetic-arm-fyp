"""Input surface overlaps, not a clinical attachment or motion validation."""
from repair_inputs import *

def solid(path):
    v,f=read_obj(str(path));v,f=weld_stl(v[f]);return mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))

def run():
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));bones={r['id']:solid(ROOT/r['source']) for r in inv['models'] if r['kind']=='bone'};rows=[]
    for r in inv['models']:
        if r['kind']!='muscle':continue
        path=OUT/'repaired'/(r['id']+'.obj')
        if not path.exists():raise RuntimeError('Missing closed muscle '+r['id'])
        s=solid(path)
        if s.status()!=mf.Error.NoError:raise RuntimeError(r['id']+' invalid repaired mesh')
        overlaps=[]
        for bid,b in bones.items():
            vol=max(0.,float((s^b).volume()));gap=float(s.min_gap(b,10));overlaps.append({'bone':bid,'intersection_mm3':vol,'gap_mm_capped_10':gap})
        rows.append({'id':r['id'],'name':r['name'],'volume_mm3':float(s.volume()),'bone_pairs':overlaps,'muscle_to_bone_transform':'identity, preserved source coordinates'})
    report={'scope':'original anatomical surfaces, detached zero-thickness triangle islands removed; no vertex displacement; not engineered parts or muscle attachments','all_sources_unchanged':all(hashlib.sha256((ROOT/r['source']).read_bytes()).hexdigest()==r['sha256'] for r in inv['models']),'muscles':rows,'interpretation':'Intersections are recorded, not automatically relieved or accepted as legal dynamic contact. Nearest bone is not proof of anatomical origin/insertion.'}
    (OUT/'input_fit_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print('source preserved',report['all_sources_unchanged'],'pairs',sum(len(x['bone_pairs']) for x in rows),'intersecting pairs',sum(p['intersection_mm3']>.001 for x in rows for p in x['bone_pairs']))

if __name__=='__main__':run()
