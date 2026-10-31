"""Actual neutral M4 manufactured meshes. Does not animate assumed TPU strain."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,os,sys,xml.etree.ElementTree as ET
import numpy as np,mujoco
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate'
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15';folder=BASE/variant
    data=json.loads((folder/'assembly.json').read_text(encoding='utf8'));design=json.loads((folder/'design.json').read_text(encoding='utf8'))
    root=ET.Element('mujoco',model='M4_actual_neutral_rigid_side_lever');ET.SubElement(root,'compiler',meshdir='.',angle='degree');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1000',offheight='900');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1');ET.SubElement(world,'light',pos='1 0 .5',dir='-1 0 -1',diffuse='.6 .6 .6')
    for r in data['instances']:
        key=r['key'];ET.SubElement(asset,'mesh',name=key,file=os.path.relpath(ROOT/r['file'],folder).replace('\\','/'),scale='.001 .001 .001')
        color='.15 .55 .8 1' if key=='side_lever' else '.85 .12 .1 1' if r['kind']=='TPU' else '.85 .82 .75 1' if key in ['FJ3362','FJ3384','FJ3368','FJ3391','FJ3290'] else '.55 .6 .65 1'
        ET.SubElement(world,'geom',name=key,type='mesh',mesh=key,pos=' '.join(str(v/1000) for v in r['translation_display_mm']),rgba=color,contype='0',conaffinity='0')
    ET.indent(root);path=folder/'neutral_view.xml';ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
    m=_load_project_model(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d);cam=mujoco.MjvCamera();detail='--detail' in sys.argv
    cam.lookat[:]=(np.array(design['joint_centre_display_mm'])+np.array(design['zero_forearm_direction'])*25)/1000 if detail else [-.05,0,.02]
    cam.distance=.38 if detail else 1.1;cam.azimuth=0 if detail else -90;cam.elevation=-10 if detail else -4
    if '--render' in sys.argv:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(m,900,1000);panel=Image.new('RGB',(2000,900),'white')
        for i,azimuth in enumerate([-30,35] if detail else [-90,90]):
            cam.azimuth=azimuth;renderer.update_scene(d,camera=cam);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1000,65),fill='white');draw.text((12,12),'M4 | actual111 neutral meshes | removable side lever',fill='black');draw.text((12,35),'No solved tendon route, TPU deformation or actuation shown',fill='black');panel.paste(im,(i*1000,0))
        output=folder/('lever_detail.png' if detail else 'whole_neutral.png');panel.save(output);renderer.close();print(output);return
    from mujoco import viewer
    import time
    with viewer.launch_passive(m,d) as v:
        v.cam.lookat[:]=cam.lookat;v.cam.distance=cam.distance;v.cam.azimuth=cam.azimuth;v.cam.elevation=cam.elevation
        while v.is_running():v.sync();time.sleep(.03)
if __name__=='__main__':run()
