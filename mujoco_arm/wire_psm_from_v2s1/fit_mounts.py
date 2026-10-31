"""Surface-rooted mounting candidates, selected by sampled steel clearance.

No TPU validity claim: offset eyes require physical brackets and muscle-volume
checks. Search does not alter bone geometry or adopt the rejected gimbal.
"""
from common import *
import itertools

def outer_triangles(name,part):
    tri=part['tri'];center=tri.mean(1);n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.maximum(np.linalg.norm(n,axis=1),1e-20)[:,None]
    if name=='scapula_mesh':
        normal=np.array([.657,-.754,0]);mask=(n@normal)>.15
    elif name in ['humerus_shaft','ulna_shaft']:
        axis=np.array([.087156,0,-.996195]);radial=center-(center@axis)[:,None]*axis
        mask=np.sum(n*radial,axis=1)>.5
    elif name=='clavicle_visual':
        axis=np.array([116.86,-94.03,-1.76]);axis/=np.linalg.norm(axis);radial=center-(center@axis)[:,None]*axis;mask=np.sum(n*radial,axis=1)>.5
    else:
        radial=center-np.array([0,0,140]);mask=np.sum(n*radial,axis=1)>.5
    return tri[mask] if mask.any() else tri

def run():
    tree=template();m,d=compile_tree(tree);parts=solids(m,d);root=tree.getroot();bybody={p['body']:(n,p) for n,p in parts.items()};routes=list(root.find('tendon'));library=pose_library();poses=[np.zeros(7)]
    for seq in library.values():poses.extend(np.array(seq)[np.unique(np.linspace(0,len(seq)-1,min(5,len(seq))).astype(int))])
    unique=np.unique(np.array(poses),axis=0);poses=np.vstack([np.zeros(7),unique[np.any(np.abs(unique)>1e-12,axis=1)]]);transforms=[]
    for q in poses:
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);transforms.append((d.xpos.copy()*1000,d.xmat.reshape(-1,3,3).copy()))
    d.qpos[:]=0;mujoco.mj_forward(m,d);mounts={};route_records=[]
    for route in routes:
        for node in route:
            if node.tag=='site':
                name=node.get('site');sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,name);body=int(m.site_bodyid[sid]);old=m.site_pos[sid]*1000
            else:
                name=node.get('geom');sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,name);body=int(m.geom_bodyid[sid]);old=m.geom_pos[sid]*1000
            bonename,part=bybody[body];foot,normal,distance=project(old,outer_triangles(bonename,part));mounts[name]={'node':node,'kind':node.tag,'id':sid,'body':body,'bone':bonename,'old':old,'foot':foot,'normal':normal}
    for rid,route in enumerate(routes):
        entries=[mounts[n.get('site') or n.get('geom')] for n in route];best=None
        # Same geometry degrees of freedom as original route: 2 ends and at
        # most one guide/passpoint, each on an explicit normal stand-off.
        trials=list(itertools.product([8.,16.,24.,32.,40.],repeat=len(entries)))
        trials.sort(key=lambda z:(max(z),sum(z)))
        for offsets in trials:
            coords=[]
            for e,offset in zip(entries,offsets):
                p=e['foot']+offset*e['normal'];coords.append(p)
                if e['kind']=='site':m.site_pos[e['id']]=p*.001
                else:
                    delta=p-m.geom_pos[e['id']]*1000;m.geom_pos[e['id']]=p*.001
                    side=e['node'].get('sidesite');sid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,side);m.site_pos[sid]+=delta*.001
            clear=0;worst=1e6;blocked=None;zero_clear=False
            for pose_index,q in enumerate(poses):
                d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);path=paths(m,d)[rid];failed=False
                for bone,part in parts.items():
                    body=part['body'];local=(path-d.xpos[body]*1000)@d.xmat[body].reshape(3,3);gap=local_gap(local,part)
                    worst=min(worst,gap)
                    if gap<.25:failed=True;blocked=(bone,float(gap));break
                if not failed:clear+=1
                if pose_index==0:
                    zero_clear=not failed
                    if failed:break
            score=(zero_clear,clear,worst,-sum(offsets))
            if best is None or score>best['score']:best={'score':score,'offsets':offsets,'coords':coords,'blocked':blocked,'coarse_clear_pose_count':clear}
            if clear==len(poses):break
        for e,p,offset in zip(entries,best['coords'],best['offsets']):
            e['selected_eye']=p;e['offset']=offset
            xml=root.find(f'.//site[@name="{e["node"].get("site")}"]') if e['kind']=='site' else root.find(f'.//geom[@name="{e["node"].get("geom")}"]')
            xml.set('pos',' '.join(str(x*.001) for x in p))
            if e['kind']=='geom':
                side=root.find(f'.//site[@name="{e["node"].get("sidesite")}"]');old=np.fromstring(side.get('pos'),sep=' ')*1000;side.set('pos',' '.join(str(x*.001) for x in old+p-e['old']))
        record={'route':route.get('name'),'coarse_pose_count':len(poses),'coarse_clear_pose_count':best['coarse_clear_pose_count'],'offsets_mm':list(best['offsets']),'first_blocked':best['blocked'],'TPU_volume_verified':False,'guide_hardware_verified':False,'anchors_frozen':False}
        route_records.append(record);print(json.dumps(record),flush=True)
        # Recompile to restore stable side sites before fitting another route.
        m,d=compile_tree(tree)
    rows=[];m,d=compile_tree(tree)
    for name,e in mounts.items():
        body=e['body'];p=e['selected_eye'];r=d.xmat[body].reshape(3,3);origin=d.xpos[body]*1000
        rows.append({'name':name,'bone':e['bone'],'body':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,body),'kind':e['kind'],'bone_root_local_mj_mm':e['foot'].tolist(),'eye_local_mj_mm':p.tolist(),'bone_root_design_mm':design(e['foot']@r.T+origin).tolist(),'eye_design_mm':design(p@r.T+origin).tolist(),'outward_normal_local_mj':e['normal'].tolist(),'standoff_mm':e['offset'],'bracket_not_yet_verified':True})
    OUT.mkdir(exist_ok=True,parents=True);tree.write(OUT/'mount_candidate.xml',encoding='utf-8');(OUT/'mount_candidate.json').write_text(json.dumps({'scope':'surface roots plus normal-offset eyes; coarse steel-only screening; not final PSM actuators','routes':route_records,'mounts':rows},indent=2))
if __name__=='__main__':run()
