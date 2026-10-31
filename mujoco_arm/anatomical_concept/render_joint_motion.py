"""Render only collision-free prefixes of the current source/SC candidate."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys
import numpy as np,mujoco
from PIL import Image,ImageDraw
from joint_motion_data import accepted_paths
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept/joint_model/SC_connection'
if '--AC' in sys.argv:OUT=OUT/('AC_compact' if '--compact' in sys.argv else 'AC_connection')
if '--GH' in sys.argv:OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept/joint_model/SC_connection/AC_compact'/('GH_preserve_head' if '--preserve-head' in sys.argv else 'GH_connection')
if '--elbow' in sys.argv:OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'

def run():
    m=_load_project_model(str(OUT/'arm.xml'));d=mujoco.MjData(m);renderer=mujoco.Renderer(m,800,800);cam=mujoco.MjvCamera();cam.lookat[:]=[0,0,.12];cam.distance=1.8;cam.azimuth=-75;cam.elevation=-5;options=mujoco.MjvOption();options.geomgroup[1]=0;options.geomgroup[2]=1;options.geomgroup[3]=1;paths=accepted_paths(OUT);frames=[]
    for name in ['elbow_flexion','GH_flexion','GH_abduction','coordinated_raise','SC_elevation']:
        seq=paths[name]
        for i in np.unique(np.linspace(0,len(seq)-1,36).astype(int)):
            d.qpos[:]=seq[i];mujoco.mj_forward(m,d);renderer.update_scene(d,camera=cam,scene_option=options);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,800,58],fill='white');draw.text((12,12),name+' | accepted sampled prefix ONLY',fill='black');draw.text((12,32),('SC/AC/GH/elbow solid candidates' if '--elbow' in sys.argv else 'SC + AC + GH candidates; elbow ideal' if '--GH' in sys.argv else 'SC + AC candidates; GH/elbow ideal' if '--AC' in sys.argv else 'SC candidate; other joints ideal')+'; NO muscle drive claim',fill='black');frames.append(im)
    frames[0].save(OUT/'joint_review_motion.gif',save_all=True,append_images=frames[1:],duration=90,loop=0);renderer.close();print('review frames',len(frames))

if __name__=='__main__':run()
