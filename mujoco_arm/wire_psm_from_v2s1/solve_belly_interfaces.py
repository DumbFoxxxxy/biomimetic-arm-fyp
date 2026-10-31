"""Solve the moving PSM end C from a fixed passive tendon length.

L_passive(C,G,B,q)=constant, C=A+s*(G-A)/|G-A|. This is an explicit
axially guided interface model, not a calibrated TPU actuator simulation.
"""
from common import *
from scipy.optimize import brentq,linprog
import csv
def run():
    spec=json.loads((OUT/'psm_module_candidate.json').read_text());modules=spec['modules'];tree=ET.parse(OUT/'psm_module_candidate.xml');root=tree.getroot();m,d=compile_tree(tree)
    ids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,z['id']+'_C') for z in modules]);A=np.array([z['belly_fixed_eye_local_mj_mm'] for z in modules]);G=np.array([z['guide_center_local_mj_mm'] for z in modules]);U=(G-A)/np.linalg.norm(G-A,axis=1)[:,None];upper=np.linalg.norm(G-A,axis=1)-24;lower=np.full(16,20.)
    qs=[];labels=[];weights=[];slices={}
    for action,seq in pose_library().items():
        start=len(qs);qs.extend(seq);labels.extend([action]*len(seq));weights.extend([1/len(seq)]*len(seq));slices[action]=[start,len(qs)]
    qs=np.array(qs);weights=np.array(weights);lo=[];hi=[]
    for q in qs:
        d.qpos[:]=np.radians(q)
        m.site_pos[ids]=(A+U*lower[:,None])*.001;mujoco.mj_forward(m,d);lo.append(d.ten_length.copy()*1000)
        m.site_pos[ids]=(A+U*upper[:,None])*.001;mujoco.mj_forward(m,d);hi.append(d.ten_length.copy()*1000)
    lo=np.array(lo);hi=np.array(hi);passive=[]
    for i in range(16):
        shared_lo=float(hi[:,i].max());shared_hi=float(lo[:,i].min())
        if shared_lo+.2<shared_hi:chosen=shared_lo+.2
        else:
            candidates=np.linspace(hi[0,i]+.01,lo[0,i]-.01,401);scores=np.array([weights@((hi[:,i]<=x)&(lo[:,i]>=x)) for x in candidates]);chosen=float(candidates[np.argmax(scores)])
        passive.append(chosen)
    passive=np.array(passive);S=np.zeros((len(qs),16));feasible=(hi<=passive)&(lo>=passive);derivative=np.zeros((len(qs),16))
    for k,q in enumerate(qs):
        d.qpos[:]=np.radians(q)
        for i,sid in enumerate(ids):
            def residual(s):
                m.site_pos[sid]=(A[i]+U[i]*s)*.001;mujoco.mj_forward(m,d);return d.ten_length[i]*1000-passive[i]
            if feasible[k,i]:s=brentq(residual,lower[i],upper[i],xtol=1e-5)
            else:s=upper[i] if hi[k,i]>passive[i] else lower[i]
            S[k,i]=s;v=residual(s);derivative[k,i]=(residual(s+.01)-v)/.01;residual(s)
        if k%100==0:print('belly interfaces',k,'/',len(qs),flush=True)
    np.savez_compressed(OUT/'belly_interfaces.npz',q_deg=qs,s_mm=S,feasible=feasible,ds_derivative=derivative,passive_lengths_mm=passive,axis=U,A=A,G=G)
    # Zero-pose geometric pull-direction check: no tension/load limits.
    d.qpos[:]=0;m.site_pos[ids]=(A+U*S[0,:,None])*.001;mujoco.mj_forward(m,d);J=np.zeros((16,m.nv))
    for j in range(m.nv):
        d.qpos[:]=0;d.qpos[j]=1e-6;mujoco.mj_forward(m,d);plus=d.ten_length.copy()
        d.qpos[j]=-1e-6;mujoco.mj_forward(m,d);J[:,j]=(plus-d.ten_length)/(2e-6)*1000
    d.qpos[:]=0;mujoco.mj_forward(m,d);moment=J/derivative[0,:,None];cone=[]
    for j in range(m.nv):
        for sign in [-1,1]:
            target=np.zeros(m.nv);target[j]=sign;fit=linprog(np.zeros(16),A_eq=moment.T,b_eq=target,bounds=(0,None),method='highs');cone.append({'joint':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j),'sign':sign,'geometric_pull_direction_exists':bool(fit.success)})
    rows=[]
    for i,z in enumerate(modules):
        valid=S[feasible[:,i],i];min_s=float(valid.min()) if len(valid) else None;max_s=float(valid.max()) if len(valid) else None;cells=max(1,int(np.ceil((max_s-20)/15))) if max_s is not None else None;nominal=20+15*cells if cells else None
        row={**z,'passive_tendon_fixed_length_mm':float(passive[i]),'C_zero_local_mj_mm':(A[i]+U[i]*S[0,i]).tolist(),'belly_zero_length_mm':float(S[0,i]),'belly_min_length_mm':min_s,'belly_max_length_mm':max_s,'reference_cell_count_15mm':cells,'reference_unloaded_length_mm':nominal,'max_cell_shortening_fraction_required':float((nominal-min_s)/(nominal-20)) if nominal else None,'kinematic_feasible_pose_count':int(feasible[:,i].sum()),'pose_count':len(qs),'belly_pose_solved':True,'actual_TPU_capability_verified':False}
        rows.append(row);root.find(f'.//site[@name="{z["id"]}_C"]').set('pos',' '.join(str(x*.001) for x in A[i]+U[i]*S[0,i]))
    tree.write(OUT/'psm_interfaces_zero.xml',encoding='utf-8')
    report={'scope':'fixed A, axis-guided C, rounded R8 guide G, bone insertion B; exact passive route-length constraint solved; physical transverse guidance and TPU properties not yet verified','modules':rows,'action_slices':slices,'zero_pose_geometric_rank':int(np.linalg.matrix_rank(moment)),'zero_pose_direction_cone':cone,'zero_pose_effective_moment_arms_mm':moment.tolist(),'joint_names':[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j) for j in range(m.njnt)],'TPU_envelope_scenario_mm':[14,34],'envelope_scenario_source':'25 percent compression rod assay only; no outer pattern or physical rating','frozen':False}
    (OUT/'psm_interfaces.json').write_text(json.dumps(report,indent=2))
    with (OUT/'moment_arm_zero.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['module',*report['joint_names']]);w.writerows([[z['id'],*moment[i]] for i,z in enumerate(modules)])
    print('solved; rank',report['zero_pose_geometric_rank'],'directions',sum(x['geometric_pull_direction_exists'] for x in cone),'/14',flush=True)
if __name__=='__main__':run()
