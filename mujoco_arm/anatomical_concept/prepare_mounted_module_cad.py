"""Prepare exact translated rigid module parts, excluding both flexible bodies."""
from build_elbow_muscle_mounts import ROOT,DEST
import json,hashlib
import numpy as np

def run():
    data=json.loads((DEST/'assembly.json').read_text(encoding='utf8'))
    folder=DEST/'rigid_module_parts';folder.mkdir(exist_ok=True)
    parts=[];sources={}
    for r in data['instances']:
        if r.get('generated') or not r.get('source_module') or r['kind']=='TPU':continue
        fid=r['source_module'];name=r['key'][len(fid)+1:]
        base=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles'/fid
        p=next(p for p in json.loads((base/'cad_input.json').read_text(encoding='utf8'))['parts'] if p['name']==name)
        p=dict(p,name=r['key'],vertices=(np.array(p['vertices'])+np.array(r['translation_display_mm'])).tolist())
        parts.append(p);sources[r['key']]={'source_OBJ':r['file'],'OBJ_sha256':hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest(),'translation_display_mm':r['translation_display_mm']}
    assert len(parts)==30
    (folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':parts}),encoding='utf8')
    (folder/'source_instances.json').write_text(json.dumps(sources,indent=2),encoding='utf8')
    print('Prepared',len(parts),'rigid module components; flexible bodies excluded')

if __name__=='__main__':run()
