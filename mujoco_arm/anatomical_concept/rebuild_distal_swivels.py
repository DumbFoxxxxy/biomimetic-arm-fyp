"""Rebuild M5B endpoints, local checks and honest rejected-direct-route evidence."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[2];SCRIPTS=ROOT/'mujoco_arm/anatomical_concept'
PY='D:/MujocoEnv/Scripts/python.exe';CAD='D:/Python3.10.6/python.exe';VIS='C:/Python314/python.exe'
def run():
    mode=[] if '--direct-reference' in sys.argv else ['--distal-load']
    commands=[(PY,'build_distal_swivels.py',[]),(CAD,'export_distal_swivels.py',[]),(PY,'check_distal_swivels.py',[]),(PY,'prepare_distal_swivel_scene.py',[]),(VIS,'view_distal_swivels.py',['--render']),(PY,'check_distal_swivel_exports.py',[]),(PY,'prepare_swivel_prints.py',[]),(CAD,'export_swivel_arm_cad.py',[]),(PY,'check_swivel_first_branch.py',[]),(PY,'screen_swivel_tendon_routes.py',[]),(PY,'check_swivel_output_direction.py',[]),(VIS,'view_distal_swivels.py',['--animate','--cutaway']),(VIS,'inspect_swivel_animation.py',[]),(PY,'write_swivel_bench_protocol.py',[]),(PY,'review_distal_swivels.py',[])]
    for python,script,args in commands:
        print('RUN',script,*args,*mode,flush=True);subprocess.run([python,str(SCRIPTS/script),*args,*mode],cwd=ROOT,check=True)
if __name__=='__main__':run()
