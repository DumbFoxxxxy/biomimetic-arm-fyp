"""Move two existing fairleads continuously to remove a dense-path wire crossing."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import minimize,linprog
from gimbal_arm_model import build_gimbal
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from shoulder_model import EXPORT,SOURCE
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances

def main():
    path=EXPORT/'full_arm_shoulder_routes_candidate.json';candidate=json.loads(path.read_text());tree,_=build_gimbal();m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m)
    routes=candidate['routes'];names=[r['id'] for r in routes];ids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,'eye_'+r) for r in names];fixed=np.array([r['fixed_eye_world_mm'] for r in routes]);from shoulder_model import v8
    target=[];seen=set()
    for seq in target_actions().values():
        for a,b in zip(seq[:-1],seq[1:]):
            for f in np.linspace(0,1,5,endpoint=False):
                z=a*(1-f)+b*f;key=tuple(np.round(z[:5],9))
                if key not in seen:seen.add(key);target.append(z)
        z=seq[-1];key=tuple(np.round(z[:5],9))
        if key not in seen:seen.add(key);target.append(z)
    ends=[];jac=[];rotations=[];angular=[];scapid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'scapula')
    for z in target:
        apply_target(m,d,z);ends.append(d.site_xpos[ids]*1000@v8.DESIGN_TO_MJ);jacs=[]
        for sid in ids:
            J=np.zeros((3,m.nv));mujoco.mj_jacSite(m,d,J,None,sid);jacs.append(v8.DESIGN_TO_MJ.T@J[:,:5])
        jac.append(jacs)
        rotations.append(v8.DESIGN_TO_MJ.T@d.xmat[scapid].reshape(3,3)@v8.DESIGN_TO_MJ)
        jr=np.zeros((3,m.nv));mujoco.mj_jacBody(m,d,None,jr,scapid);angular.append(v8.DESIGN_TO_MJ.T@jr[:,:5])
    bad=json.loads((EXPORT/'gimbal_wire_pair_report.json').read_text())['worst_pair']['wires']
    ends=np.array(ends);jac=np.array(jac);rotations=np.array(rotations);angular=np.array(angular);selected=[names.index(n) for n in bad];extra=selected[-1];base=np.r_[fixed[selected].ravel(),np.zeros(3)];keyposes=np.linspace(0,len(ends)-1,15,dtype=int)
    def get(x):
        current=fixed.copy();current[selected]=x[:6].reshape(-1,3);return current
    def margins(x):
        current=get(x);moving=ends.copy();shift=np.einsum('nij,j->ni',rotations,x[6:]);moving[:,extra]+=shift
        free=moving.copy();out=[]
        for i in selected:
            for j in range(len(names)):
                if i==j or j in selected and j>i:continue
                a=np.broadcast_to(current[i],(len(ends),3));c=np.broadcast_to(current[j],(len(ends),3))
                out.append(float(segment_distances(a,free[:,i],c,free[:,j]).min())-1.7)
        for j in range(len(names)):
            if j==extra:continue
            a=np.broadcast_to(current[j],(len(ends),3))
            out.append(float(segment_distances(moving[:,extra],moving[:,extra],a,moving[:,j]).min())-2.7)
        for k in keyposes:
            direction=current-moving[k];direction/=np.linalg.norm(direction,axis=1)[:,None];jnew=jac[k].copy();jnew[extra]+=np.cross(angular[k].T,shift[k]*.001).T
            W=np.einsum('nik,ni->kn',jnew,direction);norm=np.linalg.norm(W,axis=0);W/=norm
            # Maximize the minimum normalized positive tension.
            c=np.r_[np.zeros(len(names)),-1.];A=np.column_stack([-np.eye(len(names)),np.ones(len(names))])
            lp=linprog(c,A_ub=A,b_ub=np.zeros(len(names)),A_eq=np.block([[W,np.zeros((5,1))],[np.ones((1,len(names))),np.zeros((1,1))]]),b_eq=np.r_[np.zeros(5),1.],bounds=[(0,None)]*(len(names)+1),method='highs')
            out.append((lp.x[-1] if lp.success else -1.)-.0015)
        out.extend([float(x[6:]@w),50.-float(np.linalg.norm(x[6:]))])
        return np.array(out)
    bounds=[(-100,430),(-30,400),(-120,360)]*len(selected)+[(-40,40)]*3
    best=None
    cfg,_,_=__import__('shoulder_model').inputs();u,v,w=v8.source_basis(cfg)
    for offset in [np.zeros(9),np.r_[np.zeros(6),w*20],np.r_[30.,50.,0.,-30.,-50.,0.,w*30+u*15],np.r_[-40.,-60.,20.,40.,60.,-20.,w*25-u*20]]:
        start=np.clip(base+offset,np.array(bounds)[:,0],np.array(bounds)[:,1]);r=minimize(lambda x:np.sum(((x-base)/np.r_[np.ones(6)*100,np.ones(3)*50])**2),start,method='SLSQP',bounds=bounds,constraints=[{'type':'ineq','fun':margins}],options={'maxiter':180,'ftol':1e-7})
        okay=bool(margins(r.x).min()>-1e-5);print('CONTINUOUS',okay,r.message,'margin',margins(r.x).min(),'fixed',r.x.tolist(),flush=True)
        if okay:best=r;break
    if best is None:return
    new=get(best.x)
    for i,row in enumerate(candidate['routes']):row['fixed_eye_world_mm']=new[i].tolist()
    row=candidate['routes'][extra];row['eye_local_ac_mm']=(np.array(row['eye_local_ac_mm'])+best.x[6:]).tolist();row['additional_bracket_offset_local_mm']=(np.array(row.get('additional_bracket_offset_local_mm',[0.,0.,0.]))+best.x[6:]).tolist()
    candidate.update(continuous_fairlead_refinement=True,dense_wire_pair_constraints=True,full_geometry_and_drive_pass=False)
    path.write_text(json.dumps(candidate,indent=2))
if __name__=='__main__':main()
