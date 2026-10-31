"""Reproduce the M3 experiment, preserving its failed assembly conclusion."""

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
PY='D:/MujocoEnv/Scripts/python.exe';CAD='D:/Python3.10.6/python.exe'
def run():
    commands=[(PY,'screen_elbow_rotor_space.py',[]),(PY,'check_rotor_pocket.py',[]),(CAD,'export_rotor_pocket_native.py',['--humerus-only']),(CAD,'export_rotor_pocket_native.py',['--canonical-mesh']),(PY,'import_rotor_cad_mesh.py',[]),(CAD,'export_rotor_pocket_native.py',[]),(PY,'check_rotor_pocket.py',['--neutral']),(PY,'check_rotor_pocket.py',['--exports']),(PY,'check_rotor_pocket.py',['--audit']),(PY,'check_rotor_joint_assembly.py',[]),(PY,'check_rotor_pocket.py',['--mount-assembly']),(PY,'relieve_rotor_assembly.py',[]),(PY,'review_rotor_candidate.py',[])]
    for python,script,args in commands:
        print('RUN',script,*args,flush=True);subprocess.run([python,str(SCRIPTS/script),*args],cwd=ROOT,check=True)
if __name__=='__main__':run()
