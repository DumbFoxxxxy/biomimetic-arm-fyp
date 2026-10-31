"""View actual neutral arm or a rigid-only prescribed three-hinge path.

--render-motion creates an offscreen GIF. Interactive GUI opens only when the
user runs the default/--motion command. No TPU mesh is deformed or animated.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,sys,time
import numpy as np,mujoco
ROOT=Path(__file__).resolve().parents[2];DEST=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support'
if '--keepers' in sys.argv:DEST=DEST/'captive_keepers'

def run():
    motion='--motion' in sys.argv or '--render-motion' in sys.argv
    motion_folder=DEST.parent if '--keepers' in sys.argv else DEST
    m=_load_project_model(str(DEST/('rigid_hinges.xml' if motion else 'neutral_view.xml')));d=mujoco.MjData(m);poses=json.loads((motion_folder/'conditional_motion.json').read_text())['poses'];angles=np.array([p['elbow_deg'] for p in poses]);keys=['elbow','FJ1486_swivel','FJ1477_swivel'];values=np.array([[p['elbow_deg'],p['relative_cap_angles_deg']['FJ1486'],p['relative_cap_angles_deg']['FJ1477']] for p in poses])
    def set_angles(v):
        for key,value in zip(keys,v):jid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,key);d.qpos[m.jnt_qposadr[jid]]=np.radians(value)
        mujoco.mj_forward(m,d)
    def camera(cam):cam.lookat[:]=[-.07,0,-.08];cam.distance=.5;cam.azimuth=125;cam.elevation=-7
    if '--render-motion' in sys.argv:
        from PIL import Image,ImageDraw
        m.vis.global_.offwidth=1000;m.vis.global_.offheight=760;renderer=mujoco.Renderer(m,height=760,width=1000);cam=mujoco.MjvCamera();camera(cam);frames=[];indices=list(range(0,len(poses),6))
        if indices[-1]!=len(poses)-1:indices.append(len(poses)-1)
        indices+=indices[-2:0:-1]
        for i in indices:
            set_angles(values[i]);renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1000,62),fill='white');draw.text((10,10),f'PRESCRIBED RIGID PATH | elbow {angles[i]:.2f}deg | two cap hinges',fill='black');draw.text((10,34),'TPU omitted; no passive equilibrium / steel routing / drive simulation',fill='black');frames.append(im)
        frames[0].save(DEST/'rigid_motion.gif',save_all=True,append_images=frames[1:],duration=70,loop=0);renderer.close();print('rigid-only GIF frames',len(frames),flush=True);return
    if '--check' in sys.argv:
        mujoco.mj_forward(m,d);print('scene compiled,geoms',m.ngeom,'hingeDOF',m.nq);return
    from mujoco import viewer
    start=time.monotonic()
    with viewer.launch_passive(m,d) as v:
        camera(v.cam)
        while v.is_running():
            if motion:
                q=75*(1-np.cos((time.monotonic()-start)*.5));set_angles([np.interp(q,angles,values[:,i]) for i in range(3)])
            else:mujoco.mj_forward(m,d)
            v.sync();time.sleep(.03)

if __name__=='__main__':run()
