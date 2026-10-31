"""Fast new-yoke vs all non-forearm hard bodies, before full rigid suite."""
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation

def run():
    a=json.loads((DEST/'assembly.json').read_text());world={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in a['instances'] if r['kind']=='rigid' and r['owner'] not in ['distal_unvalidated','soft_unvalidated','elbow']};j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);s=load(DEST/'direct_muscle_yoke.obj')[0]
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in world:world['reference_'+r['id']]=load(OUT/'repaired_supports'/(r['id']+'.obj'),raw=True)[0]
    first=None;poses=[]
    for angle in range(151):
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();moved=s.transform(np.c_[R,E-R@E]);hits=[]
        for key,b in world.items():
            v=intersection(moved,b)
            if v>.001:hits.append({'part':key,'volume_mm3':v})
        if hits and first is None:first={'angle_deg':angle,'hits':hits}
        poses.append({'elbow_angle_deg':angle,'hits':hits})
    (DEST/'yoke_flexion_probe.json').write_text(json.dumps({'scope':'new direct yoke rigid geometry only, not a complete arm or TPU deformation pass','step_deg':1,'poses':poses,'first_limit':first,'frozen':False},indent=2),encoding='utf8');print('direct yoke flexion first limit',first,flush=True)

if __name__=='__main__':run()
