"""Reviewable coordinate proposal; never represents drilling authorization."""
from common import *
import csv,shutil,hashlib
from scipy.spatial.transform import Rotation

def run():
    spec=json.loads((OUT/'swivel_interfaces.json').read_text());z=np.load(OUT/'swivel_interfaces.npz');audit=json.loads((OUT/'interfaces_geometry_audit.json').read_text());m,d=compile_tree(ET.parse(OUT/'psm_module_candidate.xml'));d.qpos[:]=0;mujoco.mj_forward(m,d)
    dest=OUT/'PSM-W1';dest.mkdir(exist_ok=True);markers=[];pathsrows=[];summary=[];modules=spec['modules'];scene=[]
    for i,module in enumerate(modules):
        cb=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,module['carrier']);tb=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,module['insertion_body']);rc=d.xmat[cb].reshape(3,3);rt=d.xmat[tb].reshape(3,3);pc=d.xpos[cb]*1000;pt=d.xpos[tb]*1000
        A=np.array(module['belly_fixed_eye_local_mj_mm'])@rc.T+pc;G=np.array(module['guide_center_local_mj_mm'])@rc.T+pc;B=np.array(module['insertion_eye_local_mj_mm'])@rt.T+pt;C=np.array(module['C_zero_local_mj_mm'])@rc.T+pc;r=swivel_route(C,G,B)
        roles=[('A_root',module['carrier'],np.array(module['belly_root_local_mj_mm'])@rc.T+pc),('A_fixed',module['carrier'],A),('C0_MOBILE_NOT_HOLE',module['carrier'],C),('G_root',module['carrier'],np.array(module['guide_root_local_mj_mm'])@rc.T+pc),('G_ENTRY_PIVOT',module['carrier'],G),('O0_SHEAVE_CENTER_MOVES',module['carrier'],r['center']),('B_root',module['insertion_body'],np.array(module['insertion_root_local_mj_mm'])@rt.T+pt),('B_insertion',module['insertion_body'],B)]
        for role,owner,point in roles:markers.append([module['id'],role,owner,*design(point).tolist()])
        for kind,curve in [('active_wire',np.array([A,C])),('passive_tendon',r['path'])]:
            for number,point in enumerate(curve):pathsrows.append([module['id'],kind,number,*design(point).tolist()])
        u=(G-A)/np.linalg.norm(G-A);n=np.array(module['belly_normal_local_mj'])@rc.T;n-=u*(n@u);n/=np.linalg.norm(n);v=np.cross(n,u)
        scene.append({'id':module['id'],'A_mj_mm':A.tolist(),'C_mj_mm':C.tolist(),'G_mj_mm':G.tolist(),'B_mj_mm':B.tolist(),'belly_basis':np.column_stack([u,v,n]).tolist(),'passive_path_mj_mm':r['path'].tolist(),'guide_center_mj_mm':r['center'].tolist(),'guide_axis':r['axis'].tolist(),'guide_basis':np.column_stack([u,np.cross(r['axis'],u),r['axis']]).tolist()})
        summary.append([module['id'],module['carrier'],module['insertion_body'],module['available_AG_mm'],module['belly_zero_length_mm'],module['reference_unloaded_length_mm'],module['reference_cell_count_15mm'],module['passive_tendon_fixed_length_mm'],r['D_mm'],module['target_moment_arm_mm'],module['kinematic_feasible_pose_count'],module['scenario_feasible_pose_count']])
    for name,header,rows in [('anchors_design_mm.csv',['module','role','owner','x_mm','y_mm','z_mm'],markers),('wire_paths_zero_design_mm.csv',['module','kind','index','x_mm','y_mm','z_mm'],pathsrows),('muscle_installation.csv',['module','carrier','insertion_body','AG_install_space_mm','belly_zero_mm','reference_unloaded_mm','reference_cells','passive_tendon_mm','G_to_B_zero_mm','target_moment_arm_mm','kinematic_pose_count','scenario_pose_count'],summary)]:
        with (dest/name).open('w',encoding='utf-8',newline='') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)
    shutil.copy2(Path(__file__).with_name('interface_topology.svg'),dest/'interface_topology.svg')
    (dest/'zero_scene.json').write_text(json.dumps(scene))
    for src,target in [('swivel_interfaces.json','interfaces.json'),('interfaces_geometry_audit.json','geometry_audit.json'),('swivel_interfaces.npz','interface_motion_data.npz')]:shutil.copy2(OUT/src,dest/target)
    # Pose-range lower bounds along the existing paths; no global-max claim.
    qs=z['q_deg'];hb=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'humerus');ub=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'ulna');sb=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'scapula');v0=d.xpos[ub]-d.xpos[hb];v0/=np.linalg.norm(v0)
    w=np.cross(np.array([44.55,-106.69,112.87])-np.array([72.76,-82.21,268.45]),np.array([159.06,-7.3,187.54])-np.array([72.76,-82.21,268.45]));w/=np.linalg.norm(w)
    if w[1]>0:w=-w
    tangent=np.cross([0,0,1],w);p0=d.xpos[sb]*1000+np.array([-94.68,-103.36,-23.14])
    def metric(action,q):
        if action=='K1':return abs(q[6])
        if action=='K2':return abs(q[3])
        if action=='K3':return abs(q[4])
        if action in ['K4','K5']:return abs(q[5])
        if action.startswith('SCAP'):return abs(q[2])
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d)
        if action=='K6':quat=Rotation.from_matrix(d.xmat[sb].reshape(3,3)).as_quat();return float(np.degrees(2*np.arctan2(quat[:3]@w,quat[3])))
        if action.startswith('SLIDE'):
            p=d.xpos[sb]*1000+d.xmat[sb].reshape(3,3)@np.array([-94.68,-103.36,-23.14]);return float((p-p0)@tangent)
        v=d.xpos[ub]-d.xpos[hb];v/=np.linalg.norm(v);return float(np.degrees(np.arccos(np.clip(v0@v,-1,1))))
    ranges=[]
    for row in audit['actions']:
        action=row['action'];start,end=spec['action_slices'][action];hard=end-start
        for field in ['first_hard_geometry_failure','first_wire_pair_failure']:
            if row[field] is not None:hard=min(hard,row[field]['index'])
        h=max(0,hard-1);sc=np.flatnonzero(~np.all(z['scenario_feasible'][start:end],axis=1));operating=min(hard,int(sc[0]) if len(sc) else end-start)
        if row['first_TPU_envelope_risk'] is not None:operating=min(operating,row['first_TPU_envelope_risk']['index'])
        op=max(0,operating-1);ranges.append({'action':action,'bone_path_endpoint_metric':metric(action,qs[end-1]),'last_steel_and_C_travel_clear_sample_metric':metric(action,qs[start+h]),'conditional_25pct_compression_extension_plus_box_screen_metric':metric(action,qs[start+op]),'hard_first_failure_sample':None if hard==end-start else hard,'scenario_first_failure_sample':None if operating==end-start else operating,'unit':'mm' if action.startswith('SLIDE') else 'deg','not_global_maximum':True,'scenario_not_physical_rating':True})
    (dest/'range_screening.json').write_text(json.dumps(ranges,indent=2))
    with (dest/'range_screening.csv').open('w',encoding='utf-8-sig',newline='') as f:
        ww=csv.writer(f);ww.writerow(['action','bone_endpoint','steel_and_C_travel_prefix','conditional_TPU_scenario_prefix','unit']);ww.writerows([[r['action'],r['bone_path_endpoint_metric'],r['last_steel_and_C_travel_clear_sample_metric'],r['conditional_25pct_compression_extension_plus_box_screen_metric'],r['unit']] for r in ranges])
    manifest={'version':'PSM-W1','status':'reviewable routing and interface proposal; machining coordinates not frozen','frozen':False,'zero_pose_direction_count':sum(x['geometric_pull_direction_exists'] for x in spec['zero_pose_direction_cone']),'zero_pose_rank':spec['zero_pose_geometric_rank'],'zero_geometry':audit['zero'],'hardware_ready':False,'source_bones_unchanged':True,'TPU_actuation_verified':False,'external_motor_Bowden_harness_verified':False,'A_C_axis_guidance_hardware_verified':False,'files_sha256':{str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file() and p.name!='proposal_manifest.json'}}
    (dest/'proposal_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'coordinate_rows':len(markers),'wire_path_points':len(pathsrows),'zero_direction_count':manifest['zero_pose_direction_count'],'ranges':ranges}),flush=True)
if __name__=='__main__':run()

