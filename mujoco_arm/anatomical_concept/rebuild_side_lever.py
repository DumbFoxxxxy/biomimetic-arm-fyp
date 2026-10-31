"""Reproduce M4 CAD/checks and retain the failed direct-tendon conclusion."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import subprocess
ROOT=Path(__file__).resolve().parents[2];SCRIPTS=ROOT/'mujoco_arm/anatomical_concept'
PY='D:/MujocoEnv/Scripts/python.exe';CAD='D:/Python3.10.6/python.exe';VIS='C:/Python314/python.exe'
def run():
    commands=[(PY,'build_elbow_side_lever.py',[]),(PY,'search_side_lever_eyelets.py',[]),(PY,'build_elbow_side_lever.py',['--phase-search']),(PY,'screen_side_lever_tendons.py',[]),(PY,'check_elbow_side_lever.py',[]),(PY,'check_side_lever_assembly.py',[]),(PY,'prepare_side_lever_model.py',[]),(CAD,'export_side_lever_cad.py',[]),(VIS,'view_side_lever.py',['--render','--detail']),(VIS,'view_side_lever.py',['--render']),(PY,'audit_side_lever_model.py',['--exports']),(PY,'audit_side_lever_model.py',['--neutral']),(PY,'audit_side_lever_model.py',[]),(CAD,'export_side_lever_arm.py',[]),(PY,'check_side_lever_tendon_family.py',[]),(PY,'review_side_lever.py',[])]
    for python,script,args in commands:
        print('RUN',script,*args,flush=True);subprocess.run([python,str(SCRIPTS/script),*args],cwd=ROOT,check=True)
if __name__=='__main__':run()
