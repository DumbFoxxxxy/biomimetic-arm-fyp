"""Spatial antagonistic sector-spool design over a full 160-degree stroke."""

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
from shoulder_model import ROOT,EXPORT
sys.path.insert(0,str(ROOT.parent/'cable_design_v7'))
from elbow_spool_v7 import build_elbow_model
sys.stdout.reconfigure(encoding='utf-8')

BASE={'radius_mm':18.,'entry_radius_mm':6.,'entry_width_mm':5.,'entry_guide_mm':[0.,47.5,52.5],
      'external_anchor_mm':[0.,130.,32.5],'elbow_range_deg':[0.,160.],
      'anchor_angle_deg':0.,'entry_side_angle_deg':0.,'sector_side_angle_deg':45.,
      'fixed_sector_sides':True,'sector_side_angle_E01_deg':330.,'sector_side_angle_E02_deg':0.,
      'anchor_angle_E01_deg':0.,'anchor_angle_E02_deg':180.}

def build_stable_spool(parameters):
    tree=build_elbow_model(parameters)
    if parameters.get('fixed_sector_sides'):
        parent=tree.getroot().find('.//body[@name="humerus"]');child=tree.getroot().find('.//body[@name="ulna"]')
        for node in list(child):
            if node.tag=='site' and node.get('name','').startswith('sector_side_'):
                child.remove(node);parent.append(node)
    return tree

def evaluate(parameters,samples=161):
    tree=build_stable_spool(parameters);model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));data=mujoco.MjData(model)
    lengths=[];jac=[];wrap=[]
    for angle in np.linspace(0.,160.,samples):
        data.qpos[0]=np.radians(angle);mujoco.mj_forward(model,data)
        J=np.zeros((2,1))
        for i in range(2):
            adr=int(model.ten_J_rowadr[i]);nnz=int(model.ten_J_rownnz[i]);cols=model.ten_J_colind[adr:adr+nnz]
            J[i,cols]=data.ten_J[adr:adr+nnz]
        lengths.append(data.ten_length.copy());jac.append(J[:,0]);wrap.append(data.ten_wrapnum.tolist())
    jac=np.array(jac);lengths=np.array(lengths);opposite=bool(np.all(jac[:,0]*jac[:,1]<-1e-8));minimum=float(np.min(np.abs(jac))*1000)
    jump=float(np.max(np.abs(np.diff(lengths,axis=0)))*1000);expected_step=np.radians(160/(samples-1))*18.
    return {'pass':bool(opposite and minimum>=8. and jump<=expected_step*2.+.1),
            'opposite':opposite,'minimum_arm_mm':minimum,'maximum_length_step_mm':jump,
            'maximum_excursion_mm':float(np.ptp(lengths,axis=0).max()*1000),'parameters':parameters,
            'lengths_m':lengths.tolist(),'jacobian_m_per_rad':jac.tolist(),'wrap_point_counts':wrap}

def main():
    baseline=evaluate(BASE,1281)
    print('BASE',{k:v for k,v in baseline.items() if k not in ['lengths_m','jacobian_m_per_rad','wrap_point_counts','parameters']},flush=True)
    best=baseline
    if not baseline['pass']:
        # Vary the real route phase and preferred wrap side, rather than accept
        # a guide that changes torque sign partway through the requested stroke.
        rng=np.random.default_rng(301)
        for trial in range(1200):
            phase=float(rng.choice(np.arange(0.,360.,15.)));side1=float(rng.choice(np.arange(0.,360.,15.)));side2=float(rng.choice(np.arange(0.,360.,15.)))
            params={**BASE,'anchor_angle_E01_deg':phase,'anchor_angle_E02_deg':phase+180.,
                'sector_side_angle_E01_deg':side1,'sector_side_angle_E02_deg':side2}
            row=evaluate(params,41)
            if row['pass']:
                exact=evaluate(params,321)
                if exact['pass']:best=exact;print('FOUND',trial,flush=True);break
    EXPORT.mkdir(parents=True,exist_ok=True);(EXPORT/'axis_spool_design.json').write_text(json.dumps(best,indent=2))
    print(json.dumps({k:v for k,v in best.items() if k not in ['lengths_m','jacobian_m_per_rad','wrap_point_counts']},indent=2))
if __name__=='__main__':main()
