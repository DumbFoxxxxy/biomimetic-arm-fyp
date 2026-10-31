"""Locate real inter-body intersections by individual part; diagnostic only."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,itertools,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal
from audit_gimbal_assembly import local_solid,md
from audit_gimbal_drive import apply_target
from arm_targets import target_actions
from shoulder_model import EXPORT

def run():
    tree,_=build_gimbal();model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    parts=[]
    for g in range(model.ngeom):
        n=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,g);b=int(model.geom_bodyid[g])
        if not b or model.geom_group[g]==3 or n in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        solid,_=local_solid(model,data,g)
        if solid is not None:parts.append((n,b,solid))
    actions=target_actions();rows=[]
    for label,z in [('ZERO',actions['GH_FLEX'][0])]+[(a,actions[a][-1]) for a in ['SCAP_DOWN_ROT','SCAP_UP_ROT','GH_FLEX','GH_ABD','RAISE_170','ELBOW']]+[(f'RAISE_{i}',z) for i,z in enumerate(actions['RAISE_170'])]:
        apply_target(model,data,z)
        worlds=[(n,b,s.transform(np.column_stack([data.xmat[b].reshape(3,3),data.xpos[b]*1000]))) for n,b,s in parts]
        hits=[]
        for (na,a,sa),(nb,b,sb) in itertools.combinations(worlds,2):
            if a==b:continue
            aa=np.array(sa.bounding_box());bb=np.array(sb.bounding_box())
            if np.any(aa[:3]>=bb[3:]) or np.any(bb[:3]>=aa[3:]):continue
            inter=sa^sb;vol=inter.volume()
            if vol>1e-3:hits.append({'parts':[na,nb],'volume_mm3':float(vol),'bounds_mj_mm':list(inter.bounding_box())})
        rows.append({'pose':label,'hits':hits});print(label,len(hits),[(h['parts'],round(h['volume_mm3'],2)) for h in hits],flush=True)
    (EXPORT/'gimbal_part_intersections.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':run()
