"""Finite-diameter separation for all free steel paths at dense pose samples."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,sys,itertools,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal
from audit_gimbal_wires import wire_paths
from audit_gimbal_drive import apply_target
from arm_targets import target_actions
from shoulder_model import EXPORT,SOURCE,v8,inputs
from verification_snapshot import fingerprint
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances

def run():
    cfg,fit,anchors=inputs();tree,_=build_gimbal();m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);minimum=1e9;chest=1e9;worst=None;count=0;cache={};failures={}
    for action,seq in target_actions().items():
        points=[a*(1-f)+b*f for a,b in zip(seq[:-1],seq[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[seq[-1]]
        for index,z in enumerate(points):
            apply_target(m,d,z);paths=wire_paths(m,d);count+=1
            for i,j in itertools.combinations(range(m.ntendon),2):
                a=paths[i];b=paths[j]
                # AABB lower bound avoids most distant line/arc combinations.
                bound=np.maximum(0.,np.maximum(a.min(0)-b.max(0),b.min(0)-a.max(0)))
                if np.linalg.norm(bound)-1>max(minimum,.30):continue
                key=(i,j,np.round(a-a[0],5).tobytes(),np.round(b-a[0],5).tobytes())
                gap=cache.get(key)
                if gap is None:
                    gap=float(segment_distances(a[:-1,None,:],a[1:,None,:],b[None,:-1,:],b[None,1:,:]).min())-1.;cache[key]=gap
                if gap<minimum:minimum=gap;worst={'action':action,'index':index,'wires':[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,k) for k in [i,j]],'gap_mm':gap}
                if gap<.30 and ((i,j) not in failures or gap<failures[(i,j)]['gap_mm']):failures[(i,j)]={'action':action,'index':index,'wires':[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,k) for k in [i,j]],'gap_mm':gap}
            for path in paths:
                samples=np.vstack([a[None,:]+np.linspace(0,1,max(2,int(np.linalg.norm(b-a)/5)+1))[:,None]*(b-a)[None,:] for a,b in zip(path[:-1],path[1:])])
                chest=min(chest,float(v8.thorax_clearance(samples@v8.DESIGN_TO_MJ,cfg).min())-.5)
        print('PAIR',action,'minimum',minimum,'chest',chest,flush=True)
    report={'scope':'Finite 1 mm wires, including full terminals of distinct anchors; 5 pose samples per target segment; reconstructed <=2 degree wrap arcs; maximum arc chord error below .004 mm',
        'pose_count':count,'minimum_polyline_surface_pair_gap_mm':minimum,'minimum_analytic_chest_surface_gap_mm':chest,'worst_pair':worst,'interference_groups':list(failures.values()),'pass':minimum>=.30 and chest>=.25,'input_fingerprint':fingerprint(tree.getroot())}
    (EXPORT/'gimbal_wire_pair_report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='input_fingerprint'}),flush=True)
if __name__=='__main__':run()
