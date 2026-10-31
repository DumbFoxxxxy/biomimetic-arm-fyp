"""Reproduce the M3 experiment, preserving its failed assembly conclusion."""
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[2];SCRIPTS=ROOT/'mujoco_arm/anatomical_concept'
PY='D:/MujocoEnv/Scripts/python.exe';CAD='D:/Python3.10.6/python.exe'
def run():
    commands=[(PY,'screen_elbow_rotor_space.py',[]),(PY,'check_rotor_pocket.py',[]),(CAD,'export_rotor_pocket_native.py',['--humerus-only']),(CAD,'export_rotor_pocket_native.py',['--canonical-mesh']),(PY,'import_rotor_cad_mesh.py',[]),(CAD,'export_rotor_pocket_native.py',[]),(PY,'check_rotor_pocket.py',['--neutral']),(PY,'check_rotor_pocket.py',['--exports']),(PY,'check_rotor_pocket.py',['--audit']),(PY,'check_rotor_joint_assembly.py',[]),(PY,'check_rotor_pocket.py',['--mount-assembly']),(PY,'relieve_rotor_assembly.py',[]),(PY,'review_rotor_candidate.py',[])]
    for python,script,args in commands:
        print('RUN',script,*args,flush=True);subprocess.run([python,str(SCRIPTS/script),*args],cwd=ROOT,check=True)
if __name__=='__main__':run()
