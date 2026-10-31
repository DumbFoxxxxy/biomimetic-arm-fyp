"""Localize newly introduced AC socket restrictions against humerus sweep."""
from build_AC_connection import *
from scipy.spatial.transform import Rotation

def run():
    definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'))
    G=np.array(definition['centres_source_mm']['GH'])-SHIFT;E=np.array(definition['centres_source_mm']['elbow'])-SHIFT
    def solid(path,shift=False):
        v,f=read_obj(str(path))
        if shift:v,f=weld_stl(v[f])
        s=mf.Manifold(mf.Mesh64((v-SHIFT if shift else v).astype(np.float64),f.astype(np.uint64)))
        if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid solid '+str(path))
        return s
    inventory=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));sources={r['id']:ROOT/r['source'] for r in inventory['models']}
    scap=solid(DEST/'scapula_AC_socket.obj');oldscap=solid(sources['FJ3384'],True)
    hum=solid(sources['FJ3368'],True)
    added=scap-oldscap;retained=scap^oldscap;long=unit(E-G);axis=unit(np.cross(long,[-1,0,0]));rows=[]
    for angle in [70,71,80,90,110,119]:
        R=Rotation.from_rotvec(np.radians(angle)*axis).as_matrix();moving=hum.transform(np.c_[R,G-R@G])
        for name,s in [('retained_original_bone',retained),('added_socket_and_bridges',added)]:
            contact=s^moving;vol=max(0,float(contact.volume()))
            rows.append({'GH_abduction_deg':angle,'region':name,'overlap_mm3':vol,'intersection_bounds_display_mm':list(contact.bounding_box()) if vol>.001 else None})
    (DEST/'interference_localization.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
    for row in rows:print(row)

if __name__=='__main__':run()
