"""Integrate the existing GH/elbow meshes and measure arm directions in world.

This is geometric ROM verification, separate from tendon-driven tracking.
"""
import copy,json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import least_squares
from shoulder_model import build,inputs,EXPORT,ROOT,v8
sys.stdout.reconfigure(encoding='utf-8')

def build_arm():
    tree=build();root=tree.getroot();demo=ET.parse(ROOT.parent/'arm_tendon_demo.xml').getroot()
    scap=root.find('.//body[@name="scapula"]');gh=copy.deepcopy(demo.find('.//body[@name="gh_flex"]'))
    scap.append(gh)
    asset=root.find('asset')
    for node in demo.find('asset'):
        name=node.get('name')
        if name in ['humerus','ulna','humerus_mat','ulna_mat','joint_mat']:
            if asset.find(f"{node.tag}[@name='{name}']") is None:asset.append(copy.deepcopy(node))
    # Explicit range targets; separate front extension from flexion.
    for name,limit in [('gh_flex','-40 120'),('gh_abd','-20 120'),('gh_rot','-90 70'),('elbow','0 145')]:
        root.find(f'.//joint[@name="{name}"]').set('range',limit)
    ET.indent(tree,space='  ');return tree

def set_arm(model,data,shoulder,gh=(0.,0.,0.),elbow=0.):
    v8.set_joint_state(model,data,shoulder)
    for name,value in zip(['gh_flex','gh_abd','gh_rot','elbow'],[*gh,elbow]):
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
        data.qpos[model.jnt_qposadr[jid]]=np.radians(value)
    mujoco.mj_forward(model,data)

def world_direction(model,data):
    h=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'humerus');e=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'ulna')
    direction=v8.mj_to_design(data.xpos[e])-v8.mj_to_design(data.xpos[h]);return direction/np.linalg.norm(direction)

def mesh_samples(model,data,name):
    gid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,name);mid=int(model.geom_dataid[gid]);body=int(model.geom_bodyid[gid])
    va=int(model.mesh_vertadr[mid]);vn=int(model.mesh_vertnum[mid]);fa=int(model.mesh_faceadr[mid]);fn=int(model.mesh_facenum[mid])
    raw=model.mesh_vert[va:va+vn];world=raw@data.geom_xmat[gid].reshape(3,3).T+data.geom_xpos[gid]
    local=(world-data.xpos[body])@data.xmat[body].reshape(3,3);tri=local[model.mesh_face[fa:fa+fn]]
    return body,np.vstack([local,tri.mean(1),.5*(tri[:,0]+tri[:,1]),.5*(tri[:,1]+tri[:,2]),.5*(tri[:,2]+tri[:,0])])

def audit():
    cfg,fit,anchors=inputs();tree=build_arm();model=mujoco.MjModel.from_xml_string(ET.tostring(tree.getroot(),encoding='unicode'));data=mujoco.MjData(model)
    zero=np.array(fit['poses'][0]['q_rad']);up=[zero]+[np.array(p['q_rad']) for p in fit['poses'] if p.get('action')=='UP_ROT']
    set_arm(model,data,zero);neutral=world_direction(model,data)
    samples=[mesh_samples(model,data,name) for name in ['humerus_shaft','ulna_shaft']]
    def limb_clearance():
        values=[]
        for body,local in samples:
            world=local@data.xmat[body].reshape(3,3).T+data.xpos[body]
            values.append(float(v8.thorax_clearance(world@v8.DESIGN_TO_MJ*1000,cfg).min()))
        return min(values)
    motions=[('GH_FLEX',0,120.),('GH_EXT',0,-40.),('GH_ABD',1,120.),('GH_ADD',1,-20.),('GH_INT',2,70.),('GH_EXT_ROT',2,-90.),('ELBOW',3,145.)]
    rows=[];allposes=[]
    for name,axis,end in motions:
        minclear=1e9;endpoint={}
        for f in np.linspace(0.,1.,121):
            values=np.zeros(4);values[axis]=end*f;set_arm(model,data,zero,values[:3],values[3]);minclear=min(minclear,limb_clearance())
            allposes.append({'action':name,'fraction':float(f),'shoulder_rad':zero.tolist(),'gh_deg':values[:3].tolist(),'elbow_deg':float(values[3])})
        direction=world_direction(model,data);angle=np.degrees(np.arccos(np.clip(np.dot(direction,neutral),-1.,1.)))
        row={'action':name,'commanded_endpoint_deg':end,'world_humeral_angle_from_neutral_deg':float(angle),'min_limb_chest_clearance_mm':minclear,'geometric_chest_pass':minclear>=.5}
        rows.append(row);print(json.dumps(row),flush=True)
    # Measure actual combined raising; do not add Euler angles.
    for axis_name,axis in [('COMBINED_FLEX',0),('COMBINED_ABD',1)]:
        minclear=1e9
        for f in np.linspace(0.,1.,121):
            pos=f*(len(up)-1);k=min(int(pos),len(up)-2);blend=pos-k;shoulder=up[k]*(1-blend)+up[k+1]*blend
            gh=np.zeros(3);gh[axis]=120*f;set_arm(model,data,shoulder,gh);minclear=min(minclear,limb_clearance())
            allposes.append({'action':axis_name,'fraction':float(f),'shoulder_rad':shoulder.tolist(),'gh_deg':gh.tolist(),'elbow_deg':0.})
        angle=float(np.degrees(np.arccos(np.clip(np.dot(world_direction(model,data),neutral),-1.,1.))))
        rows.append({'action':axis_name,'gh_endpoint_deg':120.,'scapular_endpoint_deg':50.,'measured_world_raise_deg':angle,
            'min_limb_chest_clearance_mm':minclear,'geometric_chest_pass':minclear>=.5})
        print(json.dumps(rows[-1]),flush=True)
    report={'scope':'Integrated original humerus/ulna geometry and repaired shoulder; 121 samples per motion; no payload criteria; geometry only',
        'nq':model.nq,'nv':model.nv,'nu':model.nu,'motions':rows,
        'all_sampled_limb_chest_pass':all(r['geometric_chest_pass'] for r in rows),
        'limb_scapula_collision_verified':False,'gh_tendon_tracking_verified':False,'poses':allposes}
    EXPORT.mkdir(parents=True,exist_ok=True);tree.write(EXPORT/'arm_rom_kinematic.xml',encoding='utf-8',xml_declaration=True)
    (EXPORT/'arm_rom_report.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':audit()
