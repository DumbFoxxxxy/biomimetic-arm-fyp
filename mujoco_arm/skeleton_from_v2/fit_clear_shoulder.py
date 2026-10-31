"""Find shoulder pose paths preserving the repaired scapular plate."""
import json,sys
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from inspect_repairs import *
sys.path.insert(0,str(ROOT/'mujoco_arm/cable_design_v8'))
import v8_shoulder as v8

def run():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--conservative',action='store_true');args=parser.parse_args()
    m,d,s,b=load();out=ROOT/'mujoco_arm/exports/skeleton_from_v2';raw=np.load(out/'joint_relief_local_meshes.npz')
    for n in NAMES:
        verts,faces=np.unique(raw[n].reshape(-1,3),axis=0,return_inverse=True)
        s[n]=md.Manifold(md.Mesh(np.asarray(verts,dtype=np.float32),np.asarray(faces.reshape(-1,3),dtype=np.uint32)))
    points=np.asarray(s['scapula_mesh'].to_mesh().vert_properties)[:,:3];cfg=v8.load_config();sb=b['scapula_mesh'];q0=np.zeros(3)
    w=np.cross(np.array([44.55,-106.69,112.87])-np.array([72.76,-82.21,268.45]),np.array([159.06,-7.3,187.54])-np.array([72.76,-82.21,268.45]));w/=np.linalg.norm(w)
    if w[1]>0:w=-w
    def evaluate(q):
        d.qpos[:]=0;d.qpos[:3]=np.radians(q);mujoco.mj_forward(m,d);r=d.xmat[sb].reshape(3,3)
        p=points@r.T+d.xpos[sb]*1000;design=p[:,[0,2,1]]*np.array([1,1,-1])
        gap=v8.thorax_clearance(design,cfg)
        if args.conservative:
            ellipse=(np.sqrt((design[:,0]/140)**2+(design[:,2]/106)**2)-1)*106
            ellipse=np.where((design[:,1]>=0)&(design[:,1]<=280),ellipse,20)
            gap=np.minimum(gap,ellipse)
        quat=Rotation.from_matrix(r).as_quat();up=np.degrees(2*np.arctan2(quat[:3]@w,quat[3]));return gap,up
    rows=[]
    for target in np.arange(0,51,1):
        margin=.15+min(target*.15,2.)
        fit=minimize(lambda q:np.sum((q-q0)**2),q0,bounds=[(-10,55),(-25,25),(-35,35)],method='SLSQP',constraints=[{'type':'eq','fun':lambda q:evaluate(q)[1]-target},{'type':'ineq','fun':lambda q:evaluate(q)[0]-margin}],options={'ftol':1e-8,'maxiter':150})
        gap,up=evaluate(fit.x)
        row={'target_up_deg':int(target),'q_deg':fit.x.tolist(),'actual_up_deg':float(up),'requested_analytic_margin_mm':margin,'min_analytic_vertex_gap_mm':float(gap.min()),'optimizer_success':bool(fit.success)}
        if abs(up-target)>.05 or gap.min()<margin-.05:
            row['accepted_as_path_knot']=False;rows.append(row);print(json.dumps(row),flush=True);break
        row['accepted_as_path_knot']=True;rows.append(row);q0=fit.x.copy();print(json.dumps(row),flush=True)
    (out/('shoulder_path_fit_conservative.json' if args.conservative else 'shoulder_path_fit.json')).write_text(json.dumps({'scope':'analytic chest/vertices path search only; exact pair audit still required','knots':rows},indent=2))
if __name__=='__main__':run()
