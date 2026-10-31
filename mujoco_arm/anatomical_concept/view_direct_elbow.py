"""Interactive neutral elbow candidate, or rigid-only prescribed motion.

No TPU mesh is animated by an invented deformation law. GUI starts only when
the user runs this file; --check compiles both scenes without opening a window.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys,time
import numpy as np
import mujoco
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/yoke_access_relief_study'

def run():
    if '--check' in sys.argv:
        for name in ['neutral_view.xml','arm.xml']:
            m=_load_project_model(str(DEST/name));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
            print(name,'compiled; nq',m.nq,'geoms',m.ngeom)
        return
    rigid='--rigid-motion' in sys.argv
    m=_load_project_model(str(DEST/('arm.xml' if rigid else 'neutral_view.xml')));d=mujoco.MjData(m)
    print('Rigid-only 0–150deg prescribed geometry; no TPU/drive simulation.' if rigid else 'Neutral geometry only; TPU deformation and drive unverified.')
    from mujoco import viewer
    started=time.monotonic()
    with viewer.launch_passive(m,d) as v:
        v.cam.lookat[:]=[-.07,0,-.07];v.cam.distance=.55;v.cam.azimuth=125;v.cam.elevation=-7
        while v.is_running():
            if rigid:
                jid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'elbow')
                if jid<0:
                    ids=[j for j in range(m.njnt) if m.jnt_type[j]==mujoco.mjtJoint.mjJNT_HINGE]
                    assert len(ids)==1; jid=ids[0]
                d.qpos[m.jnt_qposadr[jid]]=np.radians(75*(1-np.cos((time.monotonic()-started)*.5)))
            mujoco.mj_forward(m,d);v.sync();time.sleep(.03)

if __name__=='__main__':run()
