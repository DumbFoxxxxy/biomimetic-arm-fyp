"""Dense repaired-geometry check, including connections and cable intersections."""
import json,sys
import numpy as np
import mujoco
from shoulder_model import compile_model,inputs,EXPORT,v8,input_fingerprints
sys.stdout.reconfigure(encoding='utf-8')

def mesh_local(model,data,gid):
    mid=int(model.geom_dataid[gid]); va=int(model.mesh_vertadr[mid]); vn=int(model.mesh_vertnum[mid]); fa=int(model.mesh_faceadr[mid]); fn=int(model.mesh_facenum[mid])
    raw=model.mesh_vert[va:va+vn]; R=data.geom_xmat[gid].reshape(3,3); body=int(model.geom_bodyid[gid]); B=data.xmat[body].reshape(3,3)
    points=((raw@R.T+data.geom_xpos[gid]-data.xpos[body])@B)
    triangles=points[model.mesh_face[fa:fa+fn]]
    samples=np.vstack([points,triangles.mean(1),.5*(triangles[:,0]+triangles[:,1]),.5*(triangles[:,1]+triangles[:,2]),.5*(triangles[:,2]+triangles[:,0])])
    return samples,triangles

def capsule_local(model,data,gid):
    radius,half=model.geom_size[gid,:2]; angles=np.linspace(0.,2*np.pi,24,endpoint=False)
    rings=np.array([[radius*np.cos(a),radius*np.sin(a),z] for z in np.linspace(-half,half,5) for a in angles])
    sphere=np.array([[radius*np.cos(p)*np.cos(a),radius*np.cos(p)*np.sin(a),h+radius*np.sin(p)] for h in [-half,half] for p in np.linspace(-np.pi/2,np.pi/2,9) for a in angles])
    raw=np.vstack([rings,sphere]); body=int(model.geom_bodyid[gid]);B=data.xmat[body].reshape(3,3);R=data.geom_xmat[gid].reshape(3,3)
    return (raw@R.T+data.geom_xpos[gid]-data.xpos[body])@B

def intersections(a,b,triangles):
    direction=b-a; e1=triangles[:,1]-triangles[:,0];e2=triangles[:,2]-triangles[:,0]
    h=np.cross(direction,e2);det=np.sum(e1*h,axis=1);valid=np.abs(det)>1e-12
    inv=np.where(valid,1/np.where(valid,det,1),0.);s=a-triangles[:,0]
    u=inv*np.sum(s*h,axis=1);q=np.cross(s,e1);v=inv*np.sum(direction*q,axis=1);t=inv*np.sum(e2*q,axis=1)
    return bool(np.any(valid&(u>=0)&(v>=0)&(u+v<=1)&(t>1e-6)&(t<1-1e-6)))

def audit():
    config,fit,anchors=inputs();model=compile_model();data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    objects=[]
    for gid in range(model.ngeom):
        name=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,gid) or ''
        if name in ['scapula_mesh','clavicle_mesh']:
            samples,tri=mesh_local(model,data,gid);objects.append((name,int(model.geom_bodyid[gid]),samples,tri))
        elif '_web' in name:
            objects.append((name,int(model.geom_bodyid[gid]),capsule_local(model,data,gid),None))
    fixed=np.array([a['fixed_eye_world_mm'] for a in anchors['routes']])*.001@v8.DESIGN_TO_MJ.T
    eye_ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"eye_{a['id']}") for a in anchors['routes']]
    minima={name:1e9 for name,_,_,_ in objects};minimum_cable=1e9;hits=[];pose_count=0
    zero=np.array(fit['poses'][0]['q_rad']);qsteps=[]
    for action in ['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']:
        sequence=[zero]+[np.array(p['q_rad']) for p in fit['poses'] if p.get('action')==action]
        qsteps.extend(np.degrees(np.linalg.norm(np.diff(sequence,axis=0),axis=1)).tolist())
        for j,(qa,qb) in enumerate(zip(sequence[:-1],sequence[1:])):
            for f in np.linspace(0.,1.,11):
                q=qa*(1-f)+qb*f;v8.set_joint_state(model,data,q);mujoco.mj_forward(model,data);pose_count+=1
                moving=data.site_xpos[eye_ids];frac=np.linspace(0.,1.,101)
                line=fixed[:,None,:]*(1-frac[None,:,None])+moving[:,None,:]*frac[None,:,None]
                minimum_cable=min(minimum_cable,float(v8.thorax_clearance(line.reshape(-1,3)@v8.DESIGN_TO_MJ*1000,config).min())-.5)
                for name,body,samples,tri in objects:
                    R=data.xmat[body].reshape(3,3);world=samples@R.T+data.xpos[body]
                    gap=float(v8.thorax_clearance(world@v8.DESIGN_TO_MJ*1000,config).min());minima[name]=min(minima[name],gap)
                    if tri is not None:
                        aa=(fixed-data.xpos[body])@R;bb=(moving-data.xpos[body])@R
                        for index,(a,b) in enumerate(zip(aa,bb)):
                            if intersections(a,b,tri):hits.append({'action':action,'segment':j,'fraction':float(f),'cable':index+1,'part':name})
    report={'scope':'1320 pose samples; triangle face/edge samples and capsule surfaces against inherited analytic thorax; centerline/triangle intersections',
        'pose_count':pose_count,'minimum_body_chest_clearances_mm':minima,'minimum_cable_radius_adjusted_chest_clearance_mm':minimum_cable,
        'cable_mesh_intersection_count':len(hits),'intersections_first_15':hits[:15],
        'max_adjacent_q_vector_step_deg':max(qsteps),'pass':min(minima.values())>=.5 and minimum_cable>=.5 and not hits,
        'cable_cable_clearance_verified':False,'exact_continuous_collision_certified':False,'input_sha256':input_fingerprints()}
    EXPORT.mkdir(parents=True,exist_ok=True);(EXPORT/'shoulder_geometry_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':audit()
