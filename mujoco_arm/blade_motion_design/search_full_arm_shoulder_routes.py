"""Re-search shoulder fairleads against the complete concentric-arm geometry."""
import json,sys,itertools,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog
from gimbal_arm_model import build_gimbal
from audit_gimbal_assembly import local_solid
from audit_gimbal_drive import apply_target
from arm_targets import target_actions
from shoulder_model import inputs,v8,EXPORT,SOURCE
from geometry_distance import segment_mesh_distance
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances

def main():
    cfg,fit,existing=inputs();tree,_=build_gimbal();m=mujoco.MjModel.from_xml_string(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    locals=[]
    for g in range(m.ngeom):
        n=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g);b=int(m.geom_bodyid[g])
        if not b or m.geom_group[g]==3 or n in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        s,p=local_solid(m,d,g)
        if s is not None:
            mesh=s.to_mesh();locals.append((n,b,np.asarray(mesh.vert_properties)[:,:3],np.asarray(mesh.tri_verts)))
    poses=[]
    for name,seq in target_actions().items():
        for index in [0,4,8,12,16,20]:poses.append((name,index,seq[index]))
    obstacles=[];ends=[];jac=[];eyes=np.array([r['eye_local_ac_mm'] for r in existing['routes']])
    for name,index,z in poses:
        apply_target(m,d,z);worlds=[]
        for n,b,v,f in locals:
            v=(v@d.xmat[b].reshape(3,3).T+d.xpos[b]*1000)@v8.DESIGN_TO_MJ;t=v[f]
            worlds.append((n,v.min(0),v.max(0),t,t.min(1),t.max(1)))
        obstacles.append(worlds);fk=v8.direct_kinematics(cfg,z[:5]);ends.append(eyes@fk['body_rotation'].T+fk['ac_mm']);cols=[]
        for k in range(5):
            zp=z[:5].copy();zm=zp.copy();zp[k]+=1e-6;zm[k]-=1e-6;fp=v8.direct_kinematics(cfg,zp);fm=v8.direct_kinematics(cfg,zm)
            cols.append(((eyes@fp['body_rotation'].T+fp['ac_mm'])-(eyes@fm['body_rotation'].T+fm['ac_mm']))/(2e-6))
        jac.append(np.array(cols).transpose(1,2,0))
    ends=np.array(ends);jac=np.array(jac);pool=[];fractions=np.linspace(0,1,31)
    # Repeated neutral poses are harmless; prioritize known crossing poses.
    order=sorted(range(len(poses)),key=lambda i:(poses[i][0]!='SCAP_UP_ROT',-poses[i][1]))
    for g in range(8):
        ys=[330.,380.] if g<2 else ([0.,45.] if g<4 else [120.,200.,270.,330.])
        candidates=[existing['routes'][g]['fixed_eye_world_mm']]+list(itertools.product([-80.,-40.,0.,40.,100.,180.,240.,300.,360.],ys,[-120.,-40.,40.,120.,200.,280.,340.]))
        if g==6:candidates=[existing['routes'][g]['fixed_eye_world_mm']]
        accepted=[]
        for value in candidates:
            a=np.array(value);bb=ends[:,g];line=a[None,None,:]*(1-fractions[None,:,None])+bb[:,None,:]*fractions[None,:,None]
            if v8.thorax_clearance(line.reshape(-1,3),cfg).min()<.8:continue
            okay=True
            for i in order:
                b=bb[i];b=b-5.*(b-a)/np.linalg.norm(b-a);lo=np.minimum(a,b)-.8;hi=np.maximum(a,b)+.8
                for n,pmin,pmax,t,tmin,tmax in obstacles[i]:
                    if np.any(lo>pmax) or np.any(hi<pmin):continue
                    hit=np.all(tmin<=hi,axis=1)&np.all(tmax>=lo,axis=1)
                    if hit.any() and segment_mesh_distance(a,b,t[hit])<.8:okay=False;break
                if not okay:break
            if not okay:continue
            direction=a[None,:]-bb;length=np.linalg.norm(direction,axis=1);direction/=length[:,None];w=np.einsum('nik,ni->nk',jac[:,g],direction);w/=np.linalg.norm(w,axis=1)[:,None]
            accepted.append({'fixed_eye_world_mm':a.tolist(),'w':w,'length_mm':length.tolist()})
        pool.append(accepted);print('POOL',g+1,len(accepted),flush=True)
    np.savez_compressed(EXPORT/'full_arm_route_pool.npz',**{f'g{g}_w':np.array([r['w'] for r in p]) for g,p in enumerate(pool)})
    (EXPORT/'full_arm_route_pool.json').write_text(json.dumps([[{k:v for k,v in r.items() if k!='w'} for r in p] for p in pool],indent=2))
    if any(not p for p in pool):print('NO_CLEAR_ROUTE_POOL',flush=True);return
    keyposes=[0,5,11,17,23,29,35,77,83];cache={};rng=np.random.default_rng(2911)
    def closure(selection,i):
        W=np.column_stack([pool[g][j]['w'][i] for g,j in enumerate(selection)])
        return linprog(np.zeros(8),A_eq=np.vstack([W,np.ones((1,8))]),b_eq=np.r_[np.zeros(5),1.],bounds=[(1e-5,None)]*8,method='highs').success
    def pair(g,j,h,k):
        key=(g,j,h,k)
        if key not in cache:
            a=np.broadcast_to(pool[g][j]['fixed_eye_world_mm'],(len(poses),3));c=np.broadcast_to(pool[h][k]['fixed_eye_world_mm'],(len(poses),3))
            cache[key]=bool(segment_distances(a,ends[:,g],c,ends[:,h]).min()>=1.6)
        return cache[key]
    best=None;scorebest=-1;winner=None
    for trial in range(12000):
        if best is not None and rng.random()<.85:
            selection=best.copy()
            for g in rng.choice(8,size=int(rng.integers(1,4)),replace=False):selection[g]=int(rng.integers(len(pool[g])))
        else:selection=[int(rng.integers(len(p))) for p in pool]
        nc=sum(closure(selection,i) for i in keyposes);npair=sum(pair(g,selection[g],h,selection[h]) for g in range(8) for h in range(g));score=nc+npair
        if score>scorebest or score==scorebest and rng.random()<.04:best=selection;scorebest=score;print('BEST',trial,nc,npair,flush=True)
        if nc==len(keyposes) and npair==28 and all(closure(selection,i) for i in range(len(poses))):winner=selection;break
    result={'pass':winner is not None,'geometric_pose_count':len(poses),'geometry_radius_allowance_mm':.8,'terminal_clamp_length_mm':5.,'pool_sizes':[len(p) for p in pool],'trials':trial+1,'routes':[]}
    for g,j in enumerate(winner or best):
        row={**existing['routes'][g],**{k:v for k,v in pool[g][j].items() if k!='w'}};result['routes'].append(row)
    (EXPORT/'full_arm_shoulder_routes_candidate.json').write_text(json.dumps(result,indent=2));print('RESULT',result['pass'],flush=True)
if __name__=='__main__':main()
