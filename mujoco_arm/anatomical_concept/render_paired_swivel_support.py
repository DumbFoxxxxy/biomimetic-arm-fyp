"""Actual nominal neutral meshes; passive poses/TPU are not animated."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,hashlib,sys,xml.etree.ElementTree as ET
import numpy as np,mujoco
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[2];DEST=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support'
if '--keepers' in sys.argv:DEST=DEST/'captive_keepers'

def run():
    a=json.loads((DEST/'assembly.json').read_text());root=ET.Element('mujoco',model='paired_swivel_neutral_only');ET.SubElement(root,'compiler',angle='radian');vis=ET.SubElement(root,'visual');ET.SubElement(vis,'global',offwidth='1400',offheight='900');assets=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 2',dir='0 1 -2',diffuse='.85 .85 .85');hashes={}
    for r in a['instances']:
        p=ROOT/r['file'];hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();key=r['key'];ET.SubElement(assets,'mesh',name=key,file=str(p),scale='.001 .001 .001');colour='.8 .81 .76 1'
        if r['kind']=='TPU':colour='.75 .1 .06 1'
        if key=='paired_bone_frame':colour='.2 .55 .72 1'
        if key.endswith('_cap_rotor'):colour='.8 .57 .15 1'
        ET.SubElement(world,'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,np.array(r['translation_display_mm'])/1000)),rgba=colour,contype='0',conaffinity='0')
    xml=DEST/'neutral_view.xml';ET.ElementTree(root).write(xml,encoding='utf-8',xml_declaration=True);m=_load_project_model(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d);renderer=mujoco.Renderer(m,height=900,width=1400)
    for name,lookat,distance in [('whole_arm',[-.05,0,.025],1.05),('paired_cap_detail',[-.07,0,-.07],.43)]:
        cam=mujoco.MjvCamera();cam.lookat[:]=lookat;cam.distance=distance;cam.azimuth=125;cam.elevation=-7;renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1400,64),fill='white');draw.text((12,10),f'{len(a["instances"])} actual meshes | two passive cap rotors + bone-rooted frame | MuJoCo',fill='black');draw.text((12,34),'NOMINAL NEUTRAL only; swivel angles/TPU equilibrium/steel friction/drive UNVERIFIED',fill='black');im.save(DEST/(name+'.png'))
    renderer.close();assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());(DEST/'render_inputs.json').write_text(json.dumps({'scope':'actual nominal neutral manufacturing meshes; no prescribed TPU deformation','mesh_count':len(a['instances']),'source_sha256':hashes,'assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()},indent=2),encoding='utf8');print('paired support actual neutral meshes',len(a['instances']),flush=True)

if __name__=='__main__':run()
