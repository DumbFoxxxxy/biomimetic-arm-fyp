"""Exact passive-tendon coupling for the selected swivelling sheave layout."""
from common import *
from scipy.optimize import linprog
import csv

def run():
    candidate=json.loads((OUT/'psm_module_candidate.json').read_text());modules=candidate['modules'];tree=ET.parse(OUT/'psm_module_candidate.xml');m,d=compile_tree(tree)
    if len(modules)!=16:raise RuntimeError('Incomplete module search')
    A=np.array([z['belly_fixed_eye_local_mj_mm'] for z in modules]);G=np.array([z['guide_center_local_mj_mm'] for z in modules]);B=np.array([z['insertion_eye_local_mj_mm'] for z in modules]);L=np.linalg.norm(G-A,axis=1);U=(G-A)/L[:,None]
    cb=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,z['carrier']) for z in modules]);tb=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,z['insertion_body']) for z in modules]);ids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,z['id']+'_C') for z in modules])
    def routing(q):
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);r=d.xmat.reshape(-1,3,3);p=d.xpos*1000;routes=[]
        for i in range(16):routes.append(swivel_route(A[i]@r[cb[i]].T+p[cb[i]],G[i]@r[cb[i]].T+p[cb[i]],B[i]@r[tb[i]].T+p[tb[i]]))
        return routes
    qs=[];weights=[];slices={}
    for action,seq in pose_library().items():
        start=len(qs);qs.extend(seq);weights.extend([1/len(seq)]*len(seq));slices[action]=[start,len(qs)]
    qs=np.array(qs);weights=np.array(weights);D=np.zeros((len(qs),16));angles=np.zeros_like(D)
    for k,q in enumerate(qs):
        for i,route in enumerate(routing(q)):
            D[k,i]=route['D_mm'] if route is not None else np.nan;angles[k,i]=route['wrap_deg'] if route is not None else 360
    lower=np.full(16,20.);upper=L-24;hi=L-upper+D;lo=L-lower+D;valid=np.isfinite(D)&(angles<=180+1e-6);passive=[];operating=[]
    for i in range(16):
        if not valid[0,i]:raise RuntimeError('Invalid zero guide geometry '+modules[i]['id'])
        # A declared planning scenario, NOT a measured TPU rating: choose
        # natural cell count and neutral C so zero is operable, rather than
        # "covering" motion by crushing a short belly nearly to its endblocks.
        best=None
        for cells in range(1,max(2,int(np.ceil((upper[i]-20)/15))+1)):
            nominal=20+15*cells;oplo=20+.75*(nominal-20);ophi=min(upper[i],20+1.25*(nominal-20))
            if ophi-oplo<3:continue
            for neutral in np.linspace(oplo+.15*(ophi-oplo),ophi-.15*(ophi-oplo),25):
                trial=neutral+D[:,i]-D[0,i];mask=valid[:,i]&(trial>=oplo)&(trial<=ophi);fractions=[]
                for action,(start,end) in slices.items():
                    failures=np.flatnonzero(~mask[start:end]);last=int(failures[0])-1 if len(failures) else end-start-1;fractions.append(max(0,last)/max(1,end-start-1))
                score=(min(fractions[:7]),sum(fractions),-abs(neutral-nominal)/(nominal-20))
                if best is None or score>best[0]:best=(score,float(neutral),cells,nominal,oplo,ophi)
        if best is None:raise RuntimeError('No neutral belly fits '+modules[i]['id'])
        _,neutral,cells,nominal,oplo,ophi=best;passive.append(float(L[i]+D[0,i]-neutral));operating.append({'cell_count':cells,'nominal_mm':nominal,'C_neutral_mm':neutral,'C_low_mm':oplo,'C_high_mm':ophi})
    passive=np.array(passive);S=L+D-passive;feasible=valid&(S>=lower-1e-6)&(S<=upper+1e-6)
    scenario=valid&(S>=np.array([x['C_low_mm'] for x in operating]))&(S<=np.array([x['C_high_mm'] for x in operating]))
    # The derivative ds/dq equals dD/dq by the constructed physical routing,
    # not by identifying arbitrary cable travel with TPU strain.
    moment=np.zeros((16,7))
    for j in range(7):
        q=np.zeros(7);q[j]=1e-4;plus=routing(q);q[j]=-1e-4;minus=routing(q)
        for i in range(16):moment[i,j]=-(plus[i]['D_mm']-minus[i]['D_mm'])/(2*np.radians(1e-4))
    cone=[]
    # Relative-route lengths are invariant under a common rigid transform.
    # Remove sub-micron/radian differentiation noise before cone feasibility.
    moment[np.abs(moment)<1e-4]=0.
    for j in range(7):
        for sign in [-1,1]:
            target=np.zeros(7);target[j]=sign;fit=linprog(np.zeros(16),A_eq=moment.T,b_eq=target,bounds=(0,None),method='highs');cone.append({'joint':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j),'sign':sign,'geometric_pull_direction_exists':bool(fit.success)})
    rows=[]
    for i,z in enumerate(modules):
        lengths=S[feasible[:,i],i];smin=float(lengths.min());smax=float(lengths.max());cells=operating[i]['cell_count'];nominal=operating[i]['nominal_mm']
        rows.append({**z,'passive_tendon_fixed_length_mm':float(passive[i]),'C_zero_local_mj_mm':(A[i]+U[i]*S[0,i]).tolist(),'belly_zero_length_mm':float(S[0,i]),'belly_min_length_mm':smin,'belly_max_length_mm':smax,'reference_cell_count_15mm':cells,'reference_unloaded_length_mm':nominal,'max_cell_shortening_fraction_required':float((nominal-smin)/(nominal-20)),'kinematic_feasible_pose_count':int(feasible[:,i].sum()),'scenario_feasible_pose_count':int(scenario[:,i].sum()),'operating_scenario':operating[i],'pose_count':len(qs),'belly_pose_solved':True,'actual_TPU_capability_verified':False})
        tree.getroot().find(f'.//site[@name="{z["id"]}_C"]').set('pos',' '.join(str(x*.001) for x in A[i]+U[i]*S[0,i]))
    np.savez_compressed(OUT/'swivel_interfaces.npz',q_deg=qs,s_mm=S,feasible=feasible,scenario_feasible=scenario,passive_lengths_mm=passive,axis=U,A=A,G=G,B=B,D_mm=D,wrap_deg=angles)
    report={'architecture':'motor active wire -> fixed PSM end A -> internal lattice -> moving end C (active-wire termination); separate fixed-length passive tendon C -> entry tangent G of R8 swivelling sheave -> insertion B. A and G fixed to one bone; C is not a drilling point.','coupling':'Lp=(|AG|-s)+D(G,B,q), hence s=|AG|+D-Lp; straight axial belly-end guidance assumed, not physically validated','modules':rows,'action_slices':slices,'zero_pose_geometric_rank':int(np.linalg.matrix_rank(moment)),'zero_pose_direction_cone':cone,'zero_pose_effective_moment_arms_mm':moment.tolist(),'joint_names':[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j) for j in range(7)],'TPU_envelope_scenario_mm':[14,34],'planning_strain_scenario':{'cell_compression_fraction':.25,'cell_extension_fraction':.25,'physically_calibrated':False,'extension_assay_available':False},'TPU_envelope_calibrated':False,'frozen':False}
    tree.write(OUT/'swivel_interfaces_zero.xml',encoding='utf-8');(OUT/'swivel_interfaces.json').write_text(json.dumps(report,indent=2))
    with (OUT/'swivel_moment_arm_zero.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['module',*report['joint_names']]);w.writerows([[z['id'],*moment[i]] for i,z in enumerate(modules)])
    print('rank',report['zero_pose_geometric_rank'],'positive directions',sum(x['geometric_pull_direction_exists'] for x in cone),'/14',flush=True)
if __name__=='__main__':run()
