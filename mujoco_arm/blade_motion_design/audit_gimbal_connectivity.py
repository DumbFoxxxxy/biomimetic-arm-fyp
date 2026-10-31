"""Find disconnected physical components within each rigid carrier."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal
from audit_gimbal_assembly import local_solid,md
from shoulder_model import EXPORT,v8

def run():
    tree,_=build_gimbal();m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);bs={}
    for g in range(m.ngeom):
        b=int(m.geom_bodyid[g]);n=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g)
        if not b or m.geom_group[g]==3 or n in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        s,p=local_solid(m,d,g)
        if s is not None:bs.setdefault(b,[]).append((n,s))
    rows=[]
    for b,pairs in bs.items():
        if len(pairs)<2:continue
        pieces=sorted(md.Manifold.batch_boolean([s for n,s in pairs],md.OpType.Add).decompose(),key=lambda s:-s.volume());components=[]
        for p in pieces:
            names=[n for n,s in pairs if (p^s).volume()>1e-3]
            vv=np.asarray(p.to_mesh().vert_properties)[:,:3]@v8.DESIGN_TO_MJ
            components.append({'volume_mm3':p.volume(),'bounds_design_mm':[vv.min(0).tolist(),vv.max(0).tolist()],'parts':names})
        row={'body':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,b),'positive_material_component_count':sum(p.volume()>1e-3 for p in pieces),'closed_cavity_shell_count':sum(p.volume()<-1e-3 for p in pieces),'components':components};rows.append(row);print(json.dumps(row),flush=True)
    (EXPORT/'gimbal_connectivity_report.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':run()
