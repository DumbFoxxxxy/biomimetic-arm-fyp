"""Classify first full-degree GH collision sample into retained vs added material."""
from build_GH_connection import *
from audit_joint_model import controls_to_q
import mujoco

def run():
    paths=json.loads((DEST/'motion_paths.json').read_text(encoding='utf8'))
    definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'))
    G=np.array(definition['centres_source_mm']['GH'])-SHIFT;A=np.array(definition['centres_source_mm']['AC'])-SHIFT
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']}
    original=load(ROOT/byid['FJ3384']['source'],True)[0];scap=load(DEST/'scapula_AC_GH_sockets.obj')[0];hum=load(DEST/'humerus_GH_ball.obj')[0]
    regions={'retained_original_scapula':scap^original,'added_engineering_geometry':scap-original}
    m=mujoco.MjModel.from_xml_path(str(DEST/'arm.xml'));d=mujoco.MjData(m);rows=[]
    for name in ['GH_flexion','GH_abduction']:
        pose=next((p for p in paths[name] if any(set(h['pair'])=={'FJ3384','FJ3368'} for h in p['hits'])),None)
        if not pose:continue
        d.qpos[:]=controls_to_q(m,np.array(pose['controls_deg']));mujoco.mj_forward(m,d)
        def transform(s,body,C):
            bid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,body);R=d.xmat[bid].reshape(3,3)
            return s.transform(np.c_[R,d.xpos[bid]*1000-R@C])
        moving=transform(hum,'GH',G);contacts=[]
        for region,s in regions.items():
            intersection=transform(s,'AC',A)^moving;vol=max(0,float(intersection.volume()))
            contacts.append({'region':region,'volume_mm3':vol,'intersection_bounds_display_mm':list(intersection.bounding_box()) if vol>.001 else None})
        rows.append({'probe':name,'scope':'first colliding full-degree sample, not refined boundary','controls_deg':pose['controls_deg'],'contacts':contacts})
    (DEST/'limit_localization.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
    for row in rows:print(row['probe'],[(c['region'],c['volume_mm3']) for c in row['contacts']])

if __name__=='__main__':run()
