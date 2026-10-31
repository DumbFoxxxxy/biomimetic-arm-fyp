"""Test circumferential GH cable routing around a humeral guide sleeve.

This is an explicit diagnostic model; force spanning alone is not acceptance.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog
from audit_arm_rom import build_arm,set_arm
from shoulder_model import inputs,EXPORT,v8,add_connector,tendon_jacobian,KP,KV
from search_gh_anchors import A,B,C,poses
sys.stdout.reconfigure(encoding='utf-8')

def build_wrap():
    routes=json.loads((EXPORT/'gh_anchor_search_geometry.json').read_text())['routes']
    tree=build_arm();root=tree.getroot();scap=root.find('.//body[@name="scapula"]');humerus=root.find('.//body[@name="humerus"]')
    gh=np.array([159.06,187.54,7.30]);ac=np.array(v8.load_config()['geometry_mm']['AC'])
    axis=v8.design_vector_to_mj(A);quat=v8.quat_from_z_axis(axis)
    ET.SubElement(humerus,'geom',{'name':'gh_circumferential_guide','type':'cylinder','pos':v8.vec(v8.design_to_mj(45*A)),
        'quat':v8.vec(quat),'size':'.0225 .06','contype':'0','conaffinity':'0','density':'0','rgba':'.55 .3 .8 .35'})
    for index,r in enumerate(routes):
        rid=r['id'];origin=np.array(r['origin_local_gh_mm']);end=np.array(r['insertion_local_gh_mm'])
        ET.SubElement(scap,'site',{'name':f'{rid}_root','pos':v8.vec(v8.design_to_mj(gh-ac+origin)),'size':'.0015'})
        ET.SubElement(humerus,'site',{'name':f'{rid}_insert','pos':v8.vec(v8.design_to_mj(end)),'size':'.0015'})
        angle=2*np.pi*index/8+((-1)**index)*np.pi/2
        side=45*A+45*(np.cos(angle)*B+np.sin(angle)*C)
        ET.SubElement(humerus,'site',{'name':f'{rid}_side','pos':v8.vec(v8.design_to_mj(side)),'size':'.001'})
        if r['support_root_local_gh_mm'] is not None:add_connector(scap,f'{rid}_support',gh-ac+np.array(r['support_root_local_gh_mm']),gh-ac+origin,2.)
        add_connector(humerus,f'{rid}_cuff_spoke',45*A+19.5*(end-45*A)/np.linalg.norm(end-45*A),end,1.5)
        spatial=ET.SubElement(root.find('tendon'),'spatial',{'name':rid,'width':'.001','rgba':'.25 .8 .95 1'})
        ET.SubElement(spatial,'site',{'site':f'{rid}_root'});ET.SubElement(spatial,'geom',{'geom':'gh_circumferential_guide','sidesite':f'{rid}_side'});ET.SubElement(spatial,'site',{'site':f'{rid}_insert'})
        ET.SubElement(root.find('actuator'),'position',{'name':f'pull_{rid}','tendon':rid,'kp':str(KP),'kv':str(KV),
            'ctrllimited':'true','ctrlrange':'0 2','forcelimited':'true','forcerange':'-200 0'})
    return tree

def audit():
    tree=build_wrap();model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));data=mujoco.MjData(model)
    cfg,fit,anchors=inputs();zero=np.array(fit['poses'][0]['q_rad']);rows=[]
    dofs=[model.jnt_dofadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)] for name in ['gh_flex','gh_abd','gh_rot']]
    for pose in poses():
        set_arm(model,data,zero,pose['q_deg']);J=tendon_jacobian(model,data)[8:16][:,dofs]
        norm=np.linalg.norm(J,axis=1);W=(J/np.maximum(norm[:,None],1e-10)).T
        lp=linprog(np.zeros(8),A_eq=np.vstack([W,np.ones((1,8))]),b_eq=[0,0,0,1],bounds=[(1e-5,None)]*8,method='highs')
        wrap_ranges=[]
        for j in range(8,16):
            start=int(data.ten_wrapadr[j]);count=int(data.ten_wrapnum[j])
            points=data.wrap_xpos.reshape(-1,3)[start:start+count]
            # Report native wrap path points in humerus local design coordinates.
            bid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'humerus');local=(points-data.xpos[bid])@data.xmat[bid].reshape(3,3)@v8.DESIGN_TO_MJ*1000
            wrap_ranges.append({'tendon':j,'point_count':count,'axial_min_mm':float(np.min(local@A)),'axial_max_mm':float(np.max(local@A))})
        rows.append({'id':pose['id'],'spanning':bool(lp.success),'rank':int(np.linalg.matrix_rank(W)),'wrap_path':wrap_ranges})
    result={'scope':'Eight native GH tendons with humeral cylinder wrap; no load capacity; diagnostic only',
        'pose_count':len(rows),'spanning_pass_count':sum(r['spanning'] for r in rows),'failure_ids':[r['id'] for r in rows if not r['spanning']],
        'full_geometry_verified':False,'rows':rows}
    tree.write(EXPORT/'gh_wrap_diagnostic.xml',encoding='utf-8',xml_declaration=True)
    (EXPORT/'gh_wrap_diagnostic.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':audit()
