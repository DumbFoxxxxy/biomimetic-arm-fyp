"""Blade-inspired no-payload shoulder prototype, preserving SC/AC centers.

Uses isolated mesh offsets and connecting webs; native CAD is not overwritten.
Eight pull-only tendons stabilize the 5-DOF shoulder without assumed contact.
"""
from pathlib import Path
import json,sys,hashlib,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent.parent
SOURCE=PROJECT/'tmp/shoulder_motion_design_20261029'
EXPORT=ROOT.parent/'exports/blade_motion_design'
sys.path.insert(0,str(ROOT.parent/'cable_design_v8'))
import v8_shoulder as v8

SCAP_OFFSET_MM=12.
CLAV_OFFSET_MM=np.array([0.,0.,-8.])
KP=10000.
KV=150.

def input_fingerprints():
    paths=[SOURCE/'paper_motion_fit_repaired.json',SOURCE/'anchor_search_hybrid_channels_clear.json',EXPORT/'scapula_hybrid_body_local.stl']
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def inputs():
    config=v8.load_config()
    fit=json.loads((SOURCE/'paper_motion_fit_repaired.json').read_text())
    anchors=json.loads((SOURCE/'anchor_search_hybrid_channels_clear.json').read_text())
    return config,fit,anchors

def add_connector(body,name,a_mm,b_mm,radius_mm=3.):
    a=v8.design_to_mj(np.asarray(a_mm)); b=v8.design_to_mj(np.asarray(b_mm))
    ET.SubElement(body,'geom',{'name':name,'type':'capsule','fromto':v8.vec(np.r_[a,b]),
        'size':str(radius_mm*.001),'rgba':'0.7 0.72 0.77 1','density':'500',
        'contype':'0','conaffinity':'0'})

def build():
    config,fit,anchors=inputs()
    tree=v8.build_shoulder_tree(config,0); root=tree.getroot()
    root.find('option').set('gravity','0 0 0')
    root.find('option').set('timestep','0.001')
    root.set('model','blade_paper_shoulder_no_contact')
    scap=root.find('.//body[@name="scapula"]'); clav=root.find('.//body[@name="clavicle_body"]')
    old_pad=scap.find('geom[@name="v8_scapula_contact_pad"]')
    if old_pad is not None:scap.remove(old_pad)
    w=v8.source_basis(config)[2]
    for name,delta in [('scapula_mesh',w*SCAP_OFFSET_MM),('clavicle_mesh',CLAV_OFFSET_MM)]:
        geom=root.find(f'.//geom[@name="{name}"]')
        pos=np.fromstring(geom.get('pos','0 0 0'),sep=' ')
        geom.set('pos',v8.vec(pos+v8.design_to_mj(delta)))
    channel_mesh=EXPORT/'scapula_hybrid_body_local.stl'
    if not channel_mesh.exists():raise RuntimeError('Generate the actual hybrid channel geometry first')
    root.find('.//asset/mesh[@name="scapula"]').set('file',str(channel_mesh))
    root.find('.//geom[@name="scapula_mesh"]').set('pos','0 0 0')
    add_connector(clav,'sc_web',np.zeros(3),CLAV_OFFSET_MM)
    ac=np.array(config['geometry_mm']['AC']); sc=np.array(config['geometry_mm']['SC'])
    add_connector(clav,'ac_clavicle_web',ac-sc,ac-sc+CLAV_OFFSET_MM)
    for name,p in [('ac',ac),('gh',np.array([159.06,187.54,7.30]))]:
        add_connector(scap,f'{name}_scapula_web',p-ac,p-ac+w*SCAP_OFFSET_MM)
    # Per-route stand-off is exported explicitly; connect the shifted plate
    # surface to the retained cable eye instead of silently moving joint centers.
    for row in anchors['routes']:
        moving=np.array(row['eye_local_ac_mm']); surface=np.array(row['mesh_point_local_ac_mm'])
        add_connector(scap,f"eye_web_{row['id']}",surface+w*SCAP_OFFSET_MM,moving,1.5)
        ET.SubElement(scap,'site',{'name':f"eye_{row['id']}",'pos':v8.vec(v8.design_to_mj(moving)),'size':'.0015'})
        ET.SubElement(root.find('worldbody'),'site',{'name':f"fixed_{row['id']}",
            'pos':v8.vec(v8.design_to_mj(np.array(row['fixed_eye_world_mm']))),'size':'.002'})
    tendon=root.find('tendon')
    if tendon is None:tendon=ET.SubElement(root,'tendon')
    actuator=root.find('actuator')
    if actuator is None:actuator=ET.SubElement(root,'actuator')
    for row in anchors['routes']:
        rid=row['id']; spatial=ET.SubElement(tendon,'spatial',{'name':rid,'width':'.001','rgba':'0.95 0.4 0.1 1'})
        ET.SubElement(spatial,'site',{'site':f'fixed_{rid}'})
        ET.SubElement(spatial,'site',{'site':f'eye_{rid}'})
        ET.SubElement(actuator,'position',{'name':f'pull_{rid}','tendon':rid,'kp':str(KP),'kv':str(KV),
            'forcelimited':'true','forcerange':'-200 0','ctrllimited':'true','ctrlrange':'0 2'})
    # Any pad geoms in the inherited diagnostic model are visual only.
    for geom in root.findall('.//geom'):
        geom.set('contype','0');geom.set('conaffinity','0')
    ET.indent(tree,space='  ')
    return tree

def compile_model():
    return mujoco.MjModel.from_xml_string(ET.tostring(build().getroot(),encoding='unicode'))

def q_from_native(model,data,reference=None):
    """Only used to report state; native ball quaternion remains authoritative."""
    return data.qpos.copy()

def tendon_jacobian(model,data):
    # MuJoCo can store ten_J sparsely; use mj_mulJacVec? Direct finite
    # differences in native tangent coordinates avoid representation assumptions.
    baseline=data.qpos.copy(); out=np.empty((model.ntendon,model.nv)); step=1e-6
    for k in range(model.nv):
        vel=np.zeros(model.nv);vel[k]=1.
        data.qpos[:]=baseline;mujoco.mj_integratePos(model,data.qpos,vel,step);mujoco.mj_forward(model,data);lp=data.ten_length.copy()
        data.qpos[:]=baseline;mujoco.mj_integratePos(model,data.qpos,vel,-step);mujoco.mj_forward(model,data);lm=data.ten_length.copy()
        out[:,k]=(lp-lm)/(2*step)
    data.qpos[:]=baseline;mujoco.mj_forward(model,data)
    return out

def equilibrium_tensions(model,data):
    J=tendon_jacobian(model,data)
    result=linprog(np.ones(model.ntendon),A_eq=J.T,b_eq=np.zeros(model.nv),
        bounds=[(.5,None)]*model.ntendon,method='highs')
    if not result.success:raise RuntimeError('No positive self-stress at target')
    return result.x,J

if __name__=='__main__':
    EXPORT.mkdir(parents=True,exist_ok=True)
    build().write(EXPORT/'shoulder_no_contact.xml',encoding='utf-8',xml_declaration=True)
    model=compile_model();print({'nq':model.nq,'nv':model.nv,'nu':model.nu,'ntendon':model.ntendon})
