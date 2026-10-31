"""Full pair checks of retained-shape candidate along declared action paths."""
import itertools,json,hashlib
import numpy as np
from scipy.optimize import least_squares,brentq
from scipy.spatial.transform import Rotation
from inspect_repairs import *
OUT=ROOT/'mujoco_arm/exports/skeleton_from_v2'

def load_candidate(cache='joint_relief_local_meshes.npz'):
    m,d,s,b=load();raw=np.load(OUT/cache)
    for n in NAMES:
        v,f=np.unique(raw[n].reshape(-1,3),axis=0,return_inverse=True)
        solid=md.Manifold(md.Mesh(np.asarray(v,dtype=np.float32),np.asarray(f.reshape(-1,3),dtype=np.uint32)))
        shells=solid.decompose()
        # Explicitly documented isolated cutting remnants, not joint carriers.
        chips=[p for p in shells if 0<p.volume()<=10]
        s[n]=solid-md.Manifold.batch_boolean(chips,md.OpType.Add) if chips else solid
    return m,d,s,b

def run():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--cache',default='joint_relief_local_meshes.npz');parser.add_argument('--fit',default='shoulder_path_fit.json');args=parser.parse_args()
    m,d,s,b=load_candidate(args.cache);sb=b['scapula_mesh'];hb=b['humerus_shaft'];ub=b['ulna_shaft']
    w=np.cross(np.array([44.55,-106.69,112.87])-np.array([72.76,-82.21,268.45]),np.array([159.06,-7.3,187.54])-np.array([72.76,-82.21,268.45]));w/=np.linalg.norm(w)
    if w[1]>0:w=-w
    def pose(q):d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d)
    pose(np.zeros(7));v0=d.xpos[ub]-d.xpos[hb];v0/=np.linalg.norm(v0);p0=d.xpos[sb]*1000+np.array([-94.68,-103.36,-23.14])
    slide_axis=np.cross(np.array([0.,0.,1.]),w);slide_axis/=np.linalg.norm(slide_axis)
    def metrics(q):
        pose(q);r=d.xmat[sb].reshape(3,3);v=d.xpos[ub]-d.xpos[hb];v/=np.linalg.norm(v);quat=Rotation.from_matrix(r).as_quat()
        p=d.xpos[sb]*1000+r@np.array([-94.68,-103.36,-23.14])
        return {'upper_arm_world_deg':float(np.degrees(np.arccos(np.clip(v0@v,-1,1)))),'scapula_up_twist_deg':float(np.degrees(2*np.arctan2(quat[:3]@w,quat[3]))),'contact_marker_slide_mm':float((p-p0)@slide_axis)}
    def check(q,gaps=False):
        pose(q);worlds={n:world(d,x,b[n]) for n,x in s.items()};hits=[];gaprows=[]
        for a,c in itertools.combinations(NAMES,2):
            aa=worlds[a];bb=worlds[c];ba=np.array(aa.bounding_box());bc=np.array(bb.bounding_box())
            if np.all(ba[:3]<bc[3:]) and np.all(bc[:3]<ba[3:]):
                hit=aa^bb;volume=hit.volume()
                if volume>.001:hits.append({'bones':[a,c],'volume_mm3':volume})
            if gaps:gaprows.append({'bones':[a,c],'gap_mm_capped_at_2':aa.min_gap(bb,2)})
        return hits,gaprows
    fitted=json.loads((OUT/args.fit).read_text())['knots'];shoulder=[np.zeros(7)]+[np.r_[k['q_deg'],np.zeros(4)] for k in fitted if k['accepted_as_path_knot']]
    actions={};targets={}
    for name,axis,end,target,floor in [('K1',6,145,145,116),('K2',3,120,120,96),('K3',4,120,120,96),('K4',5,70,70,56),('K5',5,-90,90,72)]:
        q=np.zeros(7);q[axis]=end;actions[name]=[np.zeros(7),q];targets[name]=(target,floor)
    actions['K6']=shoulder;targets['K6']=(50,36)
    # First raise GH in the outward plane, then add shoulder motion. Scaling
    # both motions together artificially truncates arm elevation when the
    # shoulder reaches its own first obstacle early.
    gh_raise=np.zeros(7);gh_raise[4]=120
    seq=[np.zeros(7),gh_raise]+[z+gh_raise for z in shoulder[1:]]
    actions['K7']=seq;targets['K7']=(170,132)
    old=json.loads((ROOT/'mujoco_arm/cable_design_v2/action_cable_matrix.json').read_text())['actions'];joints={mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,j):j for j in range(m.njnt)}
    for a in old:
        if a['id'].startswith('K'):continue
        q=np.zeros(7)
        for name,value in a['q_target_deg'].items():q[joints[name]]=value
        if a['id'].startswith('SLIDE'):
            requested=15 if a['id']=='SLIDE_FWD' else -15
            def displacement(angle):
                z=np.zeros(7);z[1]=angle;return metrics(z)['contact_marker_slide_mm']-requested
            q[1]=brentq(displacement,0,15) if requested>0 else brentq(displacement,-15,0)
        actions[a['id']]=[np.zeros(7),q];targets[a['id']]=(a.get('target_deg',a.get('target_mm')),None)
    rows=[];saved_paths={}
    for name,knots in actions.items():
        dense=[np.asarray(knots[0])]
        for a,c in zip(knots[:-1],knots[1:]):
            count=max(1,int(np.ceil(np.max(np.abs(c-a)))));dense.extend(a+(c-a)*f for f in np.linspace(0,1,count+1)[1:])
        accepted=[];first=None;lo=np.zeros(7);minimum_gaps={}
        for index,q in enumerate(dense):
            hits,current_gaps=check(q,True)
            if not hits:
                for gap in current_gaps:
                    key='/'.join(gap['bones']);minimum_gaps[key]=min(minimum_gaps.get(key,2),gap['gap_mm_capped_at_2'])
            if hits:
                hi=q.copy()
                while np.max(np.abs(hi-lo))>(.02 if name.startswith('SLIDE') else .1):
                    mid=(hi+lo)/2;hit,_=check(mid)
                    if hit:hi=mid;hits=hit
                    else:lo=mid
                first={'dense_index':index,'q_first_blocked_deg':hi.tolist(),'pairs':hits};accepted.append(lo.copy());break
            accepted.append(q.copy());lo=q.copy()
        actual=metrics(lo);_,gaps=check(lo,True)
        reverse_failures=sum(bool(check(q)[0]) for q in reversed(accepted))
        value=abs(lo[6]) if name=='K1' else abs(lo[3]) if name=='K2' else abs(lo[4]) if name=='K3' else abs(lo[5]) if name in ['K4','K5'] else actual['scapula_up_twist_deg'] if name=='K6' else actual['upper_arm_world_deg'] if name=='K7' else abs(lo[2]) if name.startswith('SCAP') else None
        row={'action':name,'target':targets[name][0],'floor':targets[name][1],'actual_range_deg':value,'actual_slide_mm':actual['contact_marker_slide_mm'] if name.startswith('SLIDE') else None,'endpoint_metrics':actual,'path_endpoint_clear':first is None,'first_collision':first,'accepted_forward_poses':len(accepted),'reverse_collision_count':reverse_failures,'end_pair_gaps':gaps,'minimum_sampled_pair_gaps_mm_capped_at_2':minimum_gaps,'q_last_clear_deg':lo.tolist()}
        rows.append(row);saved_paths[name]=[q.tolist() for q in accepted];print(json.dumps({k:v for k,v in row.items() if k!='end_pair_gaps'}),flush=True)
        # Save incremental evidence even if a later audit is interrupted.
        report={'frozen':False,'mesh_cache':args.cache,'mesh_cache_sha256':hashlib.sha256((OUT/args.cache).read_bytes()).hexdigest(),'native_cad_roundtrip_verified':False,'original_joint_centers_preserved':True,'collision_volume_tolerance_mm3':.001,'max_joint_sampling_step_deg':1,'first_collision_boundary_resolution_deg':.1,'continuous_collision_proof':False,'scope':'bone meshes; prescribed paths; no joint-pair exclusions; ideal kinematic joints; no force or TPU verification','actions':rows}
        (OUT/'retained_shape_rom.json').write_text(json.dumps(report,indent=2));(OUT/'retained_shape_paths.json').write_text(json.dumps(saved_paths))
    arrays={}
    for n,solid in s.items():
        mesh=solid.to_mesh();arrays[n]=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
    np.savez_compressed(OUT/'retained_shape_local_meshes.npz',**arrays)
if __name__=='__main__':run()
