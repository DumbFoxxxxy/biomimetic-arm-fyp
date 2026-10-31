"""Actual whole-arm M8 static MuJoCo view; no surrogate steel connections."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,hashlib,xml.etree.ElementTree as ET
import numpy as np,mujoco
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'mujoco_arm/exports/anatomical_concept';PARENT=OUT/'elbow_guides/ground_frame_swept';DEST=OUT/'forearm_receivers/rephased_anchors'

def run():
    screen=json.loads((DEST/'mounted_direction_screen.json').read_text());assert screen['first_rigid_limit'] is None and not screen['neutral_contacts'];tree=ET.parse(PARENT/'whole_arm.xml');root=tree.getroot();root.set('model','M8 receiving direction-study only; no complete cable route')
    for parent in root.iter():
        for child in list(parent):
            if child.tag in ['mesh','geom'] and child.get('name') in ['flexor_tendon_stop','extensor_tendon_stop']:parent.remove(child)
    root.find('./asset/mesh[@name="side_lever"]').set('file',str(DEST/'lever_receiving_spotfaces.obj'))
    for path in sorted(DEST.glob('*_neutral.obj')):
        key=path.stem;ET.SubElement(root.find('asset'),'mesh',name=key,file=str(path),scale='.001 .001 .001');colour='1 .45 .08 1' if 'crimp' in key else '.15 .65 .75 1' if 'yoke' in key else '.65 .7 .74 1';ET.SubElement(root.find('worldbody'),'geom',name=key,type='mesh',mesh=key,rgba=colour,contype='0',conaffinity='0')
    hashes={str(Path(mesh.get('file'))):hashlib.sha256(Path(mesh.get('file')).read_bytes()).hexdigest() for mesh in root.findall('./asset/mesh')};xml=DEST/'whole_arm_direction_study.xml';tree.write(xml,encoding='utf-8',xml_declaration=True);m=_load_project_model(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d);renderer=mujoco.Renderer(m,height=900,width=1400)
    for name,lookat,distance in [('whole_arm',[-.05,0,.025],1.05),('receiver_detail',[-.02,0,-.055],.25)]:
        cam=mujoco.MjvCamera();cam.lookat[:]=lookat;cam.distance=distance;cam.azimuth=125;cam.elevation=-7;renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1400,64),fill='white');draw.text((12,10),'M8 | actual receiving mechanisms + rephased lever | MuJoCo',fill='black');draw.text((12,34),'CONDITIONAL DIRECTION STUDY; upstream path incomplete; TPU fixed; NO drive simulation',fill='black');im.save(DEST/(name+'.png'))
    renderer.close();assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());(DEST/'render_inputs.json').write_text(json.dumps({'scope':'actual meshes static conditional orientation study, no cable kinks or drive imposed','mesh_count':len(root.findall('./asset/mesh')),'source_sha256':hashes,'screen_sha256':hashlib.sha256((DEST/'mounted_direction_screen.json').read_bytes()).hexdigest()},indent=2),encoding='utf8');print('M8 rendered',len(root.findall('./asset/mesh')),'mesh instances',flush=True)

if __name__=='__main__':run()
