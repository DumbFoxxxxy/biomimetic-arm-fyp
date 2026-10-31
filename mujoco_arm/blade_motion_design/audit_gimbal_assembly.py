"""Exact mesh booleans at pose checkpoints and dense surface/chest sampling.

Analytic wrap proxies are excluded from solid geometry; physical grooved
drums and bored yokes are included. No strength/payload condition is imposed.
"""

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
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from audit_geometry import mesh_local,capsule_local
from shoulder_model import EXPORT,SOURCE,inputs,v8
from verification_snapshot import fingerprint
sys.path.insert(0,str(SOURCE/'vendor'))
import manifold3d as md
sys.stdout.reconfigure(encoding='utf-8')

def local_solid(model,data,gid):
    body=int(model.geom_bodyid[gid]);B=data.xmat[body].reshape(3,3);G=data.geom_xmat[gid].reshape(3,3)
    R=B.T@G;pos=(data.geom_xpos[gid]-data.xpos[body])@B*1000;kind=model.geom_type[gid]
    if kind==mujoco.mjtGeom.mjGEOM_MESH:
        samples,tri=mesh_local(model,data,gid);tri=tri*1000;verts,faces=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True)
        solid=md.Manifold(md.Mesh(np.ascontiguousarray(verts,dtype=np.float32),np.ascontiguousarray(faces.reshape(-1,3),dtype=np.uint32)))
        return solid,samples*1000
    radius=float(model.geom_size[gid,0]*1000);half=float(model.geom_size[gid,1]*1000)
    if kind==mujoco.mjtGeom.mjGEOM_CAPSULE:
        solid=md.Manifold.cylinder(2*half,radius,circular_segments=24,center=True)
        for z in [-half,half]:solid=solid+md.Manifold.sphere(radius,circular_segments=24).translate([0,0,z])
        samples=capsule_local(model,data,gid)*1000
    elif kind==mujoco.mjtGeom.mjGEOM_CYLINDER:
        solid=md.Manifold.cylinder(2*half,radius,circular_segments=32,center=True)
        angle=np.linspace(0,2*np.pi,40,endpoint=False);raw=np.array([[radius*np.cos(a),radius*np.sin(a),z] for z in np.linspace(-half,half,5) for a in angle])
        samples=raw@R.T+pos
    elif kind==mujoco.mjtGeom.mjGEOM_SPHERE:
        solid=md.Manifold.sphere(radius,circular_segments=24);samples=np.asarray(solid.to_mesh().vert_properties)[:,:3]@R.T+pos
    else:return None,None
    return solid.transform(np.column_stack([R,pos])),samples

def audit():
    cfg,fit,anchors=inputs();tree,specs=build_gimbal();model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    parts=[];bodies={};errors=[]
    for gid in range(model.ngeom):
        name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,gid) or ''
        if model.geom_group[gid]==3 or name in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        body=int(model.geom_bodyid[gid])
        if body==0:continue
        solid,samples=local_solid(model,data,gid)
        if solid is None:continue
        if str(solid.status())!='Error.NoError':errors.append({'part':name,'status':str(solid.status())});continue
        parts.append((name,body,samples));bodies.setdefault(body,[]).append(solid)
    body_solids={body:md.Manifold.batch_boolean(solids,md.OpType.Add) for body,solids in bodies.items()}
    actions={}
    for name,seq in target_actions().items():
        actions[name]=np.array([a*(1-f)+b*f for a,b in zip(seq[:-1],seq[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[seq[-1]])
    minimum={name:1e9 for name,body,p in parts if mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,body)!='thorax'}
    analytic_worst={};overlaps=[];count=0;exact_count=0
    for action,seq in actions.items():
        for index,z in enumerate(seq):
            apply_target(model,data,z);count+=1
            for name,body,points in parts:
                if name not in minimum:continue
                world=points@data.xmat[body].reshape(3,3).T+data.xpos[body]*1000
                value=float(v8.thorax_clearance(world@v8.DESIGN_TO_MJ,cfg).min())
                if value<minimum[name]:minimum[name]=value;analytic_worst[name]={'action':action,'index':index,'clearance_mm':value}
            exact_count+=1;worlds={}
            for body,solid in body_solids.items():worlds[body]=solid.transform(np.column_stack([data.xmat[body].reshape(3,3),data.xpos[body]*1000]))
            for a,b in itertools.combinations(worlds,2):
                sa,sb=worlds[a],worlds[b];aa=np.array(sa.bounding_box());bb=np.array(sb.bounding_box())
                if np.any(aa[:3]>=bb[3:]) or np.any(bb[:3]>=aa[3:]):continue
                hit=md.Manifold.batch_boolean([sa,sb],md.OpType.Intersect);volume=float(hit.volume())
                if volume>1e-3:
                    row={'action':action,'index':index,'bodies':[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,k) for k in [a,b]],'overlap_volume_mm3':volume,'overlap_bounds_mj_mm':list(hit.bounding_box())}
                    overlaps.append(row)
        print('GEOMETRY',action,'overlap_count',len(overlaps),'minimum_analytic',min(minimum.values()),flush=True)
    report={'scope':'Full moving part surface/chest samples and exact body-union boolean checks at 5 samples per target segment; no load capacity',
        'target_knots':count,'exact_pose_checkpoints':exact_count,'invalid_meshes':errors,'minimum_part_chest_clearances_mm':minimum,
        'analytic_failures':[{**analytic_worst[name],'part':name} for name,value in minimum.items() if value<.5],
        'solid_overlap_count':len(overlaps),'overlaps':overlaps,'pass':not errors and min(minimum.values())>=.5 and not overlaps,
        'steel_path_clearance_verified':False,'input_fingerprint':fingerprint(tree.getroot()),'continuous_collision_proof':False}
    (EXPORT/'gimbal_assembly_report.json').write_text(json.dumps(report,indent=2))
    print('SUMMARY',json.dumps({'pass':report['pass'],'invalid_meshes':errors,'analytic_failures':report['analytic_failures'],'overlaps_first_10':overlaps[:10]}),flush=True)
if __name__=='__main__':audit()
