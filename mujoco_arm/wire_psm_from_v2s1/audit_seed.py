"""Audit V2 anchor/routing seeds against the actual frozen skeleton."""
from common import *
import csv,hashlib
def run():
    tree=template();m,d=compile_tree(tree);parts=solids(m,d);OUT.mkdir(exist_ok=True,parents=True);tree.write(OUT/'seed_model.xml',encoding='utf-8')
    root=tree.getroot();bybody={p['body']:p for p in parts.values()};anchor_rows=[]
    for route in root.find('tendon'):
        for node in route:
            if node.tag!='site':continue
            sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,node.get('site'));body=int(m.site_bodyid[sid]);p=m.site_pos[sid]*1000
            if body not in bybody:continue
            foot,norm,distance=project(p,bybody[body]['tri']);world=d.site_xpos[sid]*1000;root_world=foot@d.xmat[body].reshape(3,3).T+d.xpos[body]*1000
            anchor_rows.append({'route':route.get('name'),'site':node.get('site'),'body':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,body),'role':'pass_through' if 'pass_' in node.get('site') else 'muscle_endpoint_seed','old_local_mj_mm':p.tolist(),'old_global_design_mm':design(world).tolist(),'nearest_surface_root_local_mj_mm':foot.tolist(),'nearest_surface_root_design_mm':design(root_world).tolist(),'surface_normal_local_mj':norm.tolist(),'distance_to_current_surface_mm':distance,'projection_is_not_an_accepted_mount':True})
    (OUT/'seed_anchor_surface_audit.json').write_text(json.dumps(anchor_rows,indent=2))
    with (OUT/'seed_anchor_surface_audit.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['route','site','body','old_x','old_y','old_z','surface_x','surface_y','surface_z','distance_mm'])
        for row in anchor_rows:w.writerow([row['route'],row['site'],row['body'],*row['old_global_design_mm'],*row['nearest_surface_root_design_mm'],row['distance_to_current_surface_mm']])
    route_nodes=list(root.find('tendon'));owners=[]
    for route in route_nodes:
        sites=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,n.get('site')) for n in route if n.tag=='site'];owners.append((int(m.site_bodyid[sites[0]]),int(m.site_bodyid[sites[-1]])))
    rows={};cache={};lengths=[];count=0
    for action,seq in pose_library().items():
        for index,q in enumerate(seq):
            d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);wires=paths(m,d);lengths.append(d.ten_length.copy()*1000);count+=1
            for i,path in enumerate(wires):
                rid=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,i)
                for name,part in parts.items():
                    body=part['body'];local=(path-d.xpos[body]*1000)@d.xmat[body].reshape(3,3);local=local.copy()
                    # Only the 5mm terminal attachment zone of its own carrier
                    # is allowed; all other bones see the full steel path.
                    if body==owners[i][0]:
                        delta=local[1]-local[0];local[0]+=delta/max(np.linalg.norm(delta),1e-12)*min(5,np.linalg.norm(delta)*.4)
                    if body==owners[i][1]:
                        delta=local[-1]-local[-2];local[-1]-=delta/max(np.linalg.norm(delta),1e-12)*min(5,np.linalg.norm(delta)*.4)
                    key=(i,name,np.round(local,4).tobytes());gap=cache.get(key)
                    if gap is None:gap=local_gap(local,part);cache[key]=gap
                    if gap<.25 and ((i,name) not in rows or gap<rows[(i,name)]['gap_mm']):rows[(i,name)]={'route':rid,'bone':name,'gap_mm':gap,'action':action,'index':index}
        print('seed audit',action,'failure groups',len(rows),flush=True)
        (OUT/'seed_audit_partial.json').write_text(json.dumps(list(rows.values()),indent=2))
    L=np.asarray(lengths);spans=[]
    for i in range(m.ntendon):spans.append({'route':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,i),'wire_min_length_mm':float(L[:,i].min()),'wire_max_length_mm':float(L[:,i].max()),'wire_stroke_mm':float(np.ptp(L[:,i])),'TPU_shortening_verified':False})
    report={'scope':'original V2 route seeds on frozen V2-S1 meshes; all accepted bone poses; exact segment/triangle and containment tests; 1mm steel; wrap arcs reconstructed <=2deg including axial interpolation; terminal5mm allowed only on own carrier','pose_count':count,'route_count':m.ntendon,'bone_geometry_sha256':hashlib.sha256((BASE/'freeze_manifest.json').read_bytes()).hexdigest(),'wire_body_failures':list(rows.values()),'wire_length_ranges':spans,'wire_pair_clearance_verified':False,'guide_hardware_verified':False,'TPU_envelope_verified':False,'anchors_frozen':False}
    (OUT/'seed_audit.json').write_text(json.dumps(report,indent=2));print('complete',count,'failure groups',len(rows),flush=True)
if __name__=='__main__':run()
