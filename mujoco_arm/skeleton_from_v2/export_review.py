"""Export the checked mesh candidate for review, explicitly not native CAD."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import json,hashlib,xml.etree.ElementTree as ET
import numpy as np
from verify_joint_candidate import load_candidate,OUT,ROOT,NAMES,mujoco
import sys
sys.path.insert(0,str(ROOT/'mujoco_arm'))
from stl_io import write_binary_stl

def run():
    report=json.loads((OUT/'retained_shape_rom.json').read_text());paths=json.loads((OUT/'retained_shape_paths.json').read_text())
    if len(report['actions'])!=12:raise RuntimeError('All 12 actions must finish before review export')
    cache=OUT/report['mesh_cache']
    if hashlib.sha256(cache.read_bytes()).hexdigest()!=report['mesh_cache_sha256']:raise RuntimeError('Stale geometry audit')
    m,d,s,b=load_candidate(report['mesh_cache']);tree=ET.parse(ROOT/'mujoco_arm/cable_design_v2/arm_cable_design_v2.xml');root=tree.getroot();root.set('model','v2_skeleton_review_candidate_NOT_FROZEN')
    root.find('option').set('gravity','0 0 0');root.find('compiler').set('meshdir',str(OUT))
    for tag in ['tendon','actuator','sensor','contact','equality']:
        item=root.find(tag)
        if item is not None:root.remove(item)
    mapping={n:i for i,n in enumerate(NAMES)}
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='geom' and child.get('name') not in NAMES+['floor']:parent.remove(child)
    for n,solid in s.items():
        mesh=solid.to_mesh();tri=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
        path=OUT/(n+'_review_local.stl');write_binary_stl(path,tri*.001,'Checked V2 candidate; body-local metres; not native CAD')
        geom=root.find(f'.//geom[@name="{n}"]');asset=root.find(f'.//asset/mesh[@name="{geom.get("mesh")}"]')
        asset.set('file',str(path));geom.set('pos','0 0 0');geom.attrib.pop('quat',None)
    used={x.get('mesh') for x in root.findall('.//geom') if x.get('mesh')}
    for asset in list(root.find('asset')):
        if asset.tag=='mesh' and asset.get('name') not in used:root.find('asset').remove(asset)
    # Retain sites as seed references, hidden in review; not accepted muscle anchors.
    for site in root.findall('.//site'):site.set('rgba','0 0 0 0')
    ET.indent(tree,space='  ');tree.write(OUT/'skeleton_review.xml',encoding='utf-8')
    check=_load_project_model(str(OUT/'skeleton_review.xml'))
    manifest={'frozen':False,'native_cad_roundtrip_verified':False,'review_model_dofs':check.nv,'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'skeleton_review.xml',OUT/'retained_shape_rom.json',OUT/'retained_shape_paths.json']+list(OUT.glob('*_review_local.stl'))}}
    (OUT/'review_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
if __name__=='__main__':run()
