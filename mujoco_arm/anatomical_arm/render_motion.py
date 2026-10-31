"""Animate only accepted BONE prefixes; neutral TPU specimens hidden on motion."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json
import numpy as np,mujoco
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_arm/A1'

def run():
    m=_load_project_model(str(OUT/'arm.xml'));d=mujoco.MjData(m);m.vis.global_.offwidth=800;m.vis.global_.offheight=800
    renderer=mujoco.Renderer(m,800,800);cam=mujoco.MjvCamera();cam.lookat[:]=[.16,-.05,-.05];cam.distance=1.9;cam.azimuth=85;cam.elevation=-5
    paths=json.loads((OUT/'motion_paths.json').read_text(encoding='utf8'));audit=json.loads((OUT/'range_audit.json').read_text(encoding='utf8'));frames=[];opt=mujoco.MjvOption();opt.geomgroup[1]=0;opt.geomgroup[2]=0;opt.geomgroup[3]=0
    for aid in ['K1','K2','K3','K4','K5','K6','K7']:
        seq=paths[aid];first=next((i for i,p in enumerate(seq) if p['hits']),len(seq));seq=seq[:first]
        if not seq:continue
        seq=[seq[i] for i in np.unique(np.linspace(0,len(seq)-1,14).astype(int))]
        for p in seq+seq[-2::-1]:
            d.qpos[:]=np.radians(p['q_deg']);mujoco.mj_forward(m,d);renderer.update_scene(d,camera=cam,scene_option=opt);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,800,53],fill='white');draw.text((15,12),'A1 '+aid+' | accepted sampled BONE path | ideal joints',fill='black');draw.text((15,30),'Hardware installation and TPU deformation NOT validated',fill='black');frames.append(im)
    frames[0].save(OUT/'bone_motion.gif',save_all=True,append_images=frames[1:],duration=80,loop=0);renderer.close();print('frames',len(frames))

if __name__=='__main__':run()
