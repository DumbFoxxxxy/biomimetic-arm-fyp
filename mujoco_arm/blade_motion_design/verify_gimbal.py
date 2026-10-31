"""Run or reuse only passing audits for the exact current motion model."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import copy,json,subprocess,sys
from project_runtime import WorkspacePath as Path
from gimbal_arm_model import build_gimbal
from shoulder_model import EXPORT,ROOT
from verification_snapshot import fingerprint

def main():
    tree,_=build_gimbal();base=fingerprint(tree.getroot());gravity=copy.deepcopy(tree.getroot());gravity.find('option').set('gravity','0 0 -9.81');gf=fingerprint(gravity)
    steps=[('gimbal_assembly_report','audit_gimbal_assembly.py',[],base),('gimbal_wire_geometry_report','audit_gimbal_wires.py',[],base),('gimbal_wire_pair_report','audit_gimbal_wire_pairs.py',[],base),('gimbal_inverse_drive_report','audit_gimbal_inverse.py',[],base),('gimbal_gravity_drive_report','audit_gimbal_inverse.py',['--gravity'],gf)]
    for report_name,script,args,fp in steps:
        path=EXPORT/f'{report_name}.json';old=json.loads(path.read_text()) if path.exists() else {}
        is_weight_study=report_name=='gimbal_gravity_drive_report'
        if old.get('input_fingerprint')==fp and (old.get('pass') or is_weight_study):print('REUSE',report_name,'pass',old.get('pass'),flush=True);continue
        print('RUN',script,*args,flush=True);subprocess.run([sys.executable,str(ROOT/script),*args],check=True);result=json.loads(path.read_text())
        if result.get('input_fingerprint')!=fp or not result.get('pass') and not is_weight_study:raise RuntimeError(f'Current gate failed: {report_name}')
    for script in ['audit_gimbal_connectivity.py','export_gimbal_design.py']:
        print('RUN',script,flush=True);subprocess.run([sys.executable,str(ROOT/script)],check=True)
    subprocess.run([sys.executable,str(ROOT/'run_gimbal_demo.py'),'--headless','--seconds','26'],check=True)
    print('IDEAL_TENDON_REFERENCE_GATES_PASS; TPU_MUSCLE_ACTUATION_NOT_MODELED',flush=True)
if __name__=='__main__':main()
