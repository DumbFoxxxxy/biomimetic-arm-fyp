"""Reproduce M4 CAD/checks and retain the failed direct-tendon conclusion."""
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[2];SCRIPTS=ROOT/'mujoco_arm/anatomical_concept'
PY='D:/MujocoEnv/Scripts/python.exe';CAD='D:/Python3.10.6/python.exe';VIS='C:/Python314/python.exe'
def run():
    commands=[(PY,'build_elbow_side_lever.py',[]),(PY,'search_side_lever_eyelets.py',[]),(PY,'build_elbow_side_lever.py',['--phase-search']),(PY,'screen_side_lever_tendons.py',[]),(PY,'check_elbow_side_lever.py',[]),(PY,'check_side_lever_assembly.py',[]),(PY,'prepare_side_lever_model.py',[]),(CAD,'export_side_lever_cad.py',[]),(VIS,'view_side_lever.py',['--render','--detail']),(VIS,'view_side_lever.py',['--render']),(PY,'audit_side_lever_model.py',['--exports']),(PY,'audit_side_lever_model.py',['--neutral']),(PY,'audit_side_lever_model.py',[]),(CAD,'export_side_lever_arm.py',[]),(PY,'check_side_lever_tendon_family.py',[]),(PY,'review_side_lever.py',[])]
    for python,script,args in commands:
        print('RUN',script,*args,flush=True);subprocess.run([python,str(SCRIPTS/script),*args],cwd=ROOT,check=True)
if __name__=='__main__':run()
