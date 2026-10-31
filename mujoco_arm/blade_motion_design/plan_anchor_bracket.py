"""Plan a connected S09 bracket around all swept free steel paths."""
import itertools,json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal
from audit_gimbal_wires import wire_paths
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from shoulder_model import EXPORT,SOURCE,inputs,v8
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances

def main():
    cfg,fit,_=inputs();routes=json.loads((EXPORT/'full_arm_shoulder_routes_candidate.json').read_text())['routes'];row=next(r for r in routes if r['id']=='S09');end=np.array(row['eye_local_ac_mm']);start=end-np.array(row['additional_bracket_offset_local_mm']);lo=np.minimum(start,end)-35;hi=np.maximum(start,end)+35
    tree,_=build_gimbal();m=mujoco.MjModel.from_xml_string(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'scapula');segments=[];rots=[];origins=[]
    for seq in target_actions().values():
        for z in [a*(1-f)+b*f for a,b in zip(seq[:-1],seq[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[seq[-1]]:
            apply_target(m,d,z);R=d.xmat[sid].reshape(3,3);rots.append(v8.DESIGN_TO_MJ.T@R@v8.DESIGN_TO_MJ);origins.append(v8.mj_to_design(d.xpos[sid]))
            for tid,p in enumerate(wire_paths(m,d)):
                rid=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,tid)
                if rid in ['S05','S09']:
                    p=p.copy();delta=p[-1]-p[-2];p[-1]-=5*delta/np.linalg.norm(delta)
                local=(p-d.xpos[sid]*1000)@R@v8.DESIGN_TO_MJ;a=local[:-1];b=local[1:];keep=np.all(np.minimum(a,b)<=hi,axis=1)&np.all(np.maximum(a,b)>=lo,axis=1)
                if keep.any():segments.append(np.stack([a[keep],b[keep]],axis=1))
    segments=np.concatenate(segments);rot=np.array(rots);origin=np.array(origins);best=None
    for label,point in [('start',start),('end',end)]:print('ENDPOINT',label,float(segment_distances(point[None,:],point[None,:],segments[:,0],segments[:,1]).min())-2.,flush=True)
    values=sorted(itertools.product([-30.,-20.,-10.,0.,10.,20.,30.],repeat=3),key=lambda o:np.linalg.norm(o))
    for offset in values:
        delta=np.array(offset);path=np.array([start,start+delta,end+delta,end])
        if np.linalg.norm(delta)<1e-8:continue
        gap=min(float(segment_distances(a[None,:],b[None,:],segments[:,0],segments[:,1]).min())-2. for a,b in zip(path[:-1],path[1:]))
        if gap<.30:continue
        samples=np.vstack([a[None,:]+np.linspace(0,1,15)[:,None]*(b-a)[None,:] for a,b in zip(path[:-1],path[1:])]);world=np.einsum('pj,nij->npi',samples,rot)+origin[:,None,:];chest=float(v8.thorax_clearance(world.reshape(-1,3),cfg).min())-1.5
        if chest<.6:continue
        length=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum());result={'pass':True,'path_local_ac_mm':path.tolist(),'bracket_radius_mm':1.5,'minimum_swept_free_wire_surface_gap_mm':gap,'minimum_analytic_chest_gap_mm':chest,'length_mm':length,'pose_count':len(rot),'solid_assembly_audit_pass':False}
        if best is None or length<best['length_mm']:best=result;print('BRACKET',best,flush=True)
    if best is None:print('NO_BRACKET_PATH',flush=True);return
    (EXPORT/'s09_bracket_path.json').write_text(json.dumps(best,indent=2))
if __name__=='__main__':main()
