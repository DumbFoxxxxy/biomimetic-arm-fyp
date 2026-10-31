"""Read-only V2 skeletal geometry audit; no actuation or strength acceptance."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,sys,hashlib,itertools
import numpy as np
import mujoco
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'mujoco_arm/blade_motion_design'),str(ROOT/'tmp/shoulder_motion_design_20261029/vendor')]
from audit_gimbal_assembly import local_solid
OUT=ROOT/'mujoco_arm/exports/skeleton_from_v2'

def run():
    source=ROOT/'mujoco_arm/cable_design_v2';m=_load_project_model(str(source/'arm_cable_design_v2.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    selected=['thorax_geom','clavicle_visual','scapula_mesh','humerus_shaft','ulna_shaft'];parts=[]
    for name in selected:
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,name)
        if g<0:raise RuntimeError('Missing V2 bone geometry: '+name)
        s,_=local_solid(m,d,g);parts.append((name,int(m.geom_bodyid[g]),s))
    jid={mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j):j for j in range(m.njnt)}
    h=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'humerus');u=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'ulna');s=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'scapula')
    reference=d.xpos[u]-d.xpos[h];reference/=np.linalg.norm(reference);r0=d.xmat[s].reshape(3,3).copy()
    actions=json.loads((source/'action_cable_matrix.json').read_text())['actions'];actions=[{**a,'audit_id':'V2_'+a['id']} for a in actions]
    for aid,joint,end in [('K1','elbow',145.),('K2','gh_flex',120.),('K3','gh_abd',120.),('K4','gh_rot',70.),('K5','gh_rot',-90.)]:actions.append({'audit_id':'FULL_'+aid,'q_target_deg':{joint:end}})
    rows=[]
    for action in actions:
        overlaps=[];violations=[]
        for index,f in enumerate(np.linspace(0,1,21)):
            d.qpos[:]=0.
            for name,value in action['q_target_deg'].items():
                j=jid[name];d.qpos[m.jnt_qposadr[j]]=np.radians(value*f)
                if m.jnt_limited[j] and (d.qpos[m.jnt_qposadr[j]]<m.jnt_range[j,0]-1e-9 or d.qpos[m.jnt_qposadr[j]]>m.jnt_range[j,1]+1e-9):violations.append({'index':index,'joint':name,'value_deg':value*f})
            mujoco.mj_forward(m,d);world=[]
            for name,b,solid in parts:world.append((name,solid.transform(np.column_stack([d.xmat[b].reshape(3,3),d.xpos[b]*1000]))))
            for (an,aa),(bn,bb) in itertools.combinations(world,2):
                ab=np.array(aa.bounding_box());bc=np.array(bb.bounding_box())
                if np.any(ab[:3]>=bc[3:]) or np.any(bc[:3]>=ab[3:]):continue
                hit=aa^bb;volume=float(hit.volume())
                if volume>.001:overlaps.append({'index':index,'fraction':float(f),'bones':[an,bn],'volume_mm3':volume,'bounds_mj_mm':list(hit.bounding_box())})
        direction=d.xpos[u]-d.xpos[h];direction/=np.linalg.norm(direction)
        row={'action':action['audit_id'],'commanded_endpoint_deg':action['q_target_deg'],'actual_world_upper_arm_angle_from_zero_deg':float(np.degrees(np.arccos(np.clip(reference@direction,-1,1)))),'scapula_total_rotation_deg':float(np.degrees(Rotation.from_matrix(d.xmat[s].reshape(3,3)@r0.T).magnitude())),
            'sample_count':21,'joint_limit_violations':violations,'intersecting_pose_count':len(set(r['index'] for r in overlaps)),'overlaps':overlaps,'geometric_clear_at_all_samples':not overlaps and not violations}
        rows.append(row);print({k:v for k,v in row.items() if k not in ['overlaps','joint_limit_violations']},flush=True)
    OUT.mkdir(parents=True,exist_ok=True)
    report={'scope':'Unmodified V2 bone visual meshes; prescribed kinematic poses only; exact mesh booleans at 21 samples per path. Existing joint contact/intersections are reported, not automatically ignored. No forces, cable/TPU feasibility, or global maximum-ROM claim.',
        'source_xml':str(source/'arm_cable_design_v2.xml'),'source_xml_sha256':hashlib.sha256((source/'arm_cable_design_v2.xml').read_bytes()).hexdigest(),'native_velocity_dofs':m.nv,'actions':rows,'baseline_frozen':False}
    (OUT/'v2_bone_geometry_baseline.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':run()
