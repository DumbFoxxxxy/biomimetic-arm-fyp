"""Register independent access-relief geometry; preserve the previous candidate."""
from build_direct_muscle_yoke import *

def run():
    assert '--access-relief' in sys.argv and SHIFTED
    source=OUT/'direct_muscle_interfaces/mechanical_yoke_shifted'
    a=json.loads((source/'assembly.json').read_text())
    record=json.loads((DEST/'study.json').read_text())
    assert record['source_sha256']==hashlib.sha256((source/'direct_muscle_yoke.obj').read_bytes()).hexdigest()
    old=load(source/'direct_muscle_yoke.obj')[0];new=load(DEST/'direct_muscle_yoke.obj')[0]
    assert max(0,float((new-old).volume()))<.001
    assert len([s for s in new.decompose() if s.volume()>.001])==1
    for row in a['instances']:
        if row['key']=='direct_muscle_yoke':row['file']=str((DEST/'direct_muscle_yoke.obj').relative_to(ROOT))
    a.update(scope='M9 direct-muscle independent bolt-access relief candidate; soft deformation and complete assembly unverified',parent_assembly=str((source/'assembly.json').relative_to(ROOT)))
    (DEST/'assembly.json').write_text(json.dumps(a,indent=2),encoding='utf8')
    parts=json.loads((source/'cad_input.json').read_text())['parts']
    revised=record['trials'][-1]['mesh']
    parts=[revised if p['name']=='direct_muscle_yoke' else p for p in parts]
    (DEST/'cad_input.json').write_text(json.dumps({'units':'mm','parts':parts},indent=2),encoding='utf8')
    print('registered',len(a['instances']),'instances;',len(parts),'CAD parts; removal-only yoke revision',flush=True)

if __name__=='__main__':run()
