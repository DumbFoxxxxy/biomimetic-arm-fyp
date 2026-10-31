"""Real bolt-head access subtraction, retaining the original cap load plate.

The candidate is independent. A one-component gate rejects a disconnected
brace; no claim of strength or screwdriver access follows from this geometry.
"""
from build_direct_muscle_yoke import *

def run():
    source=OUT/'direct_muscle_interfaces/mechanical_yoke_shifted'
    dest=OUT/'direct_muscle_interfaces/yoke_access_relief_study';dest.mkdir(exist_ok=True)
    a=json.loads((source/'assembly.json').read_text());row=next(r for r in a['instances'] if r['key']=='direct_muscle_yoke')
    original=load(ROOT/row['file'])[0]
    d=json.loads((OUT/'direct_muscle_interfaces/shifted_muscle/FJ1486/design.json').read_text())
    F=np.array(d['basis_columns_display']);x,y,z=F.T
    placement=np.array(next(p for p in a['placements'] if p['id']=='FJ1486')['translation_display_mm'])
    pad=np.array(d['pad_centres_display_mm'][1])+placement
    q=pad+x*7
    trials=[]
    for radius in [3.5,3.8,4.0]:
        cutter=rod(q+z*6.15,q+z*42,radius,48)
        revised=(original-cutter).as_original().simplify(.001)
        count=len([s for s in revised.decompose() if s.volume()>.001])
        record={'clearance_radius_mm':radius,'cap_plate_preserved_below_local_z_mm':6.15,'components':count,'removed_mm3':float((original-revised).volume())}
        trials.append(record);print(record,flush=True)
        if count==1:
            record['mesh']=write(dest/'direct_muscle_yoke.stl',revised)
            (dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[record['mesh']]},indent=2),encoding='utf8')
            break
    (dest/'study.json').write_text(json.dumps({'source_OBJ':row['file'],'source_sha256':hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest(),'trials':trials,'strength_verified':False,'selected_as_main':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
