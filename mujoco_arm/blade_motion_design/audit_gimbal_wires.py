"""Explicit steel path audit, including reconstructed cylinder wrapping arcs."""

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
from shoulder_model import EXPORT,v8,inputs
from geometry_distance import segment_mesh_distance
from verification_snapshot import fingerprint

def wire_paths(m,d,clip_terminals=False):
    xyz=d.wrap_xpos.reshape(-1,3)*1000;obj=d.wrap_obj.ravel();out=[]
    for i in range(m.ntendon):
        a=int(d.ten_wrapadr[i]);n=int(d.ten_wrapnum[i]);p=xyz[a:a+n];o=obj[a:a+n]
        if n==2:
            if clip_terminals:
                delta=p[-1]-p[-2];p=p.copy();p[-1]-=delta/np.linalg.norm(delta)*5.
            out.append(p);continue
        p0,p1,p2,p3=p;g=int(o[1]);normal=d.geom_xmat[g].reshape(3,3)[:,2];origin=d.geom_xpos[g]*1000
        r1=p1-origin;r1-=normal*np.dot(r1,normal);r2=p2-origin;r2-=normal*np.dot(r2,normal)
        angle=np.arctan2(normal@np.cross(r1,r2),r1@r2);radius=m.geom_size[g,0]*1000
        length=d.ten_length[i]*1000-np.linalg.norm(p1-p0)-np.linalg.norm(p3-p2)
        theta=min([angle,angle+2*np.pi,angle-2*np.pi],key=lambda x:abs(abs(x)*radius-length))
        count=max(2,int(np.ceil(abs(theta)/np.radians(2.)))+1);t=np.linspace(0,theta,count)
        rr=r1[None,:]*np.cos(t)[:,None]+np.cross(normal,r1)[None,:]*np.sin(t)[:,None]
        arc=origin+normal*np.dot(p1-origin,normal)+rr
        out.append(np.vstack([p0,arc,p3]))
    return out

def audit():
    cfg,fit,anchors=inputs();tree,_=build_gimbal();m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    parts=[]
    for g in range(m.ngeom):
        name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '';body=int(m.geom_bodyid[g])
        if not body or m.geom_group[g]==3 or name in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        s,p=local_solid(m,d,g)
        if s is None:continue
        mesh=s.to_mesh();v=np.asarray(mesh.vert_properties)[:,:3];f=np.asarray(mesh.tri_verts);tri=v[f]
        parts.append((name,body,s,v.min(0),v.max(0),tri,tri.min(1),tri.max(1)))
    worst={};minimum=1e9;wire_count=0;cache={}
    for action,points in target_actions().items():
        points=np.array([a*(1-f)+b*f for a,b in zip(points[:-1],points[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[points[-1]])
        for index,z in enumerate(points):
            apply_target(m,d,z);paths=wire_paths(m,d);wire_count+=len(paths)
            for i,path in enumerate(paths):
                rid=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,i);lo=path.min(0)-.5;hi=path.max(0)+.5
                for name,body,solid,pmin,pmax,tri,tlo,thi in parts:
                    # Winding on the matching physical drum is intentional;
                    # groove planes, root radius and flange gaps are explicit.
                    if '_' in rid and name==rid.split('_')[0]+'_drum':continue
                    local=(path-d.xpos[body]*1000)@d.xmat[body].reshape(3,3)
                    if rid.startswith('S') and mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,body)=='scapula':
                        delta=local[-1]-local[-2];local=local.copy();local[-1]-=delta/np.linalg.norm(delta)*5.
                    if np.any(local.min(0)-.75>pmax) or np.any(local.max(0)+.75<pmin):continue
                    cachekey=(rid,name,np.round(local,5).tobytes());gap=cache.get(cachekey)
                    if gap is None:
                        gap=1e6
                        for a,b in zip(local[:-1],local[1:]):
                            selected=np.all(tlo<=np.maximum(a,b)+.75,axis=1)&np.all(thi>=np.minimum(a,b)-.75,axis=1)
                            if selected.any():gap=min(gap,segment_mesh_distance(a,b,tri[selected])-.5)
                            if gap<0:break
                            mid=(a+b)/2
                            if np.all(mid>pmin) and np.all(mid<pmax) and len(solid.ray_cast((pmin-np.array([123.,151.,179.])).tolist(),mid.tolist()))%2:
                                gap=-.5;break
                        cache[cachekey]=gap
                    minimum=min(minimum,gap);key=(rid,name)
                    if gap<.25 and (key not in worst or gap<worst[key]['gap_mm']):worst[key]={'wire':rid,'part':name,'gap_mm':gap,'action':action,'index':index}
        print('WIRES',action,'failures',len(worst),'min',minimum,flush=True)
        (EXPORT/'gimbal_wire_geometry_partial.json').write_text(json.dumps(list(worst.values()),indent=2))
    report={'scope':'5 samples per target segment; 1 mm steel; explicit wrap arcs at <=2 degree arc steps; final 5 mm clamp contact allowed only against the termination carrier; full wires checked against other bodies; corresponding drum winding contact excluded; radius-based triangle filters and midpoint containment checks',
        'evaluated_wire_paths':wire_count,'minimum_detected_gap_mm':minimum,'interference_groups':list(worst.values()),'pass':not worst,'pair_clearance_verified':False,'input_fingerprint':fingerprint(tree.getroot())}
    (EXPORT/'gimbal_wire_geometry_report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='input_fingerprint'}),flush=True)
if __name__=='__main__':audit()
