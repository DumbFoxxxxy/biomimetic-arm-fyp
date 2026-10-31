"""Interactive neutral anatomy/lattice review. No soft-body dynamics.
Run C:/Python314/python.exe view_arm.py; --render creates a neutral PNG.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,time,argparse
import numpy as np
import mujoco
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_arm/A1'

def run():
    parser=argparse.ArgumentParser();parser.add_argument('--render',action='store_true');parser.add_argument('--check',action='store_true');parser.add_argument('--elbow',action='store_true');parser.add_argument('--mounts',action='store_true');parser.add_argument('--fittings',action='store_true');args=parser.parse_args()
    target=OUT/'elbow_connection/muscle_mounts/end_fittings' if args.fittings else OUT/'elbow_connection/muscle_mounts' if args.mounts else OUT/'elbow_connection' if args.elbow else OUT
    m=_load_project_model(str(target/'arm.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    if args.check:print('nq',m.nq,'geoms',m.ngeom,'mesh',m.nmesh);return
    def camera(cam):cam.lookat[:]=[.15,-.05,-.25];cam.distance=1.12;cam.azimuth=85;cam.elevation=-5
    if args.render:
        from PIL import Image,ImageDraw
        m.vis.global_.offwidth=1000;m.vis.global_.offheight=960
        renderer=mujoco.Renderer(m,960,1000);cam=mujoco.MjvCamera();camera(cam);opt=mujoco.MjvOption();opt.geomgroup[2]=0;renderer.update_scene(d,camera=cam,scene_option=opt);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1000,65),fill='white');draw.text((20,15),'A1 | anatomical bones + generated TPU elbow lattices',fill='black');draw.text((20,36),'Elbow connected candidate; shoulder references remain; no soft simulation' if args.elbow else 'Neutral design candidate | blue hardware references NOT installed | no soft simulation',fill='black');im.save(target/'neutral_review.png');renderer.close();return
    from mujoco import viewer as viewer_api
    def key(k):
        if k in (ord('M'),ord('m')):viewer.opt.geomgroup[1]^=1
        if k in (ord('H'),ord('h')):viewer.opt.geomgroup[3]^=1
        if k in (ord('L'),ord('l')):viewer.opt.geomgroup[2]^=1
        if k in (ord('R'),ord('r')):d.qpos[:]=0;camera(viewer.cam)
    with viewer_api.launch_passive(m,d,key_callback=key) as viewer:
        camera(viewer.cam);viewer.opt.geomgroup[2]=0
        while viewer.is_running():
            mujoco.mj_forward(m,d)
            # Neutral rigid lattice would falsely look like a deforming muscle.
            if np.max(np.abs(d.qpos))>1e-6:viewer.opt.geomgroup[1]=0;viewer.opt.geomgroup[2]=0;viewer.opt.geomgroup[4]=0
            viewer.sync();time.sleep(.02)

if __name__=='__main__':run()
