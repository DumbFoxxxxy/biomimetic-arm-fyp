"""Expose the GH orientation singularity separately from joint actuation rank."""

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
from arm_targets import target_actions
from audit_gimbal_drive import apply_target,native_jacobian
from shoulder_model import EXPORT

if __name__=='__main__':
    tree,_=build_gimbal();m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);z=next(iter(target_actions().values()))[0].copy();hid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'humerus');columns=[int(m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n)]) for n in ['gh_flex','gh_abd','gh_rot']];rows=[]
    for angle in [0.,60.,89.,90.,91.,120.]:
        z[6]=np.radians(angle);apply_target(m,d,z);J=np.zeros((3,m.nv));mujoco.mj_jacBody(m,d,None,J,hid);singular=np.linalg.svd(J[:,columns],compute_uv=False)
        rows.append({'gh_abduction_deg':angle,'orientation_jacobian_singular_values':singular.tolist(),'orientation_rank':int(np.linalg.matrix_rank(J[:,columns],tol=1e-8)),'full_tendon_joint_jacobian_rank':int(np.linalg.matrix_rank(native_jacobian(m,d),tol=1e-8))})
    report={'scope':'The three hinges remain separately actuated; the hand orientation Jacobian loses a direction at the Euler/gimbal singularity. Passing a prescribed abduction path does not prove arbitrary independent orientation control at 90 deg.','rows':rows}
    (EXPORT/'gimbal_orientation_singularity.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
