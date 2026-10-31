"""Play full-arm tendon controls; initialize state once and integrate thereafter."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import argparse,json,time,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal
from verification_snapshot import fingerprint
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from shoulder_model import EXPORT

def main():
    actions=list(target_actions());parser=argparse.ArgumentParser();parser.add_argument('--action',choices=actions,default='RAISE_170');parser.add_argument('--seconds',type=float,default=0.);parser.add_argument('--headless',action='store_true');parser.add_argument('--gravity',action='store_true');parser.add_argument('--force-limit',type=float,default=200.);args=parser.parse_args()
    tree,_=build_gimbal()
    if args.gravity:tree.getroot().find('option').set('gravity','0 0 -9.81')
    for actuator in tree.getroot().find('actuator'):actuator.set('forcerange',f'-{args.force_limit:g} 0')
    label='gimbal_gravity' if args.gravity else 'gimbal_inverse'
    if args.force_limit!=200.:label+=f'_{args.force_limit:g}N'
    report=json.loads((EXPORT/f'{label}_drive_report.json').read_text())
    if not report['pass'] or report.get('input_fingerprint')!=fingerprint(tree.getroot()):raise RuntimeError('Drive report is failed or stale: run audit_gimbal_inverse.py for these exact inputs')
    tables=np.load(EXPORT/f'{label}_control_tables.npz');m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);apply_target(m,d,target_actions()[args.action][0]);viewer=None
    if not args.headless:
        from mujoco import viewer as mjviewer
        viewer=mjviewer.launch_passive(m,d)
    try:
        steps=0;tgrid=tables['time_s'];controls=tables[f'{args.action}_ctrl'];cycle=float(tgrid[-1])
        while (viewer is None or viewer.is_running()) and (not args.seconds or d.time<args.seconds):
            before=time.perf_counter();t=d.time%cycle;x=t/cycle*(len(tgrid)-1);i=min(int(x),len(tgrid)-2);f=x-i;d.ctrl[:]=controls[i]*(1-f)+controls[i+1]*f;mujoco.mj_step(m,d);steps+=1
            if viewer is not None:
                if steps%20==0:viewer.sync()
                time.sleep(max(0,m.opt.timestep-(time.perf_counter()-before)))
        print({'action':args.action,'seconds':float(d.time),'steps':steps,'finite':bool(np.isfinite(d.qpos).all()),'tendons':m.ntendon})
    finally:
        if viewer is not None:viewer.close()
if __name__=='__main__':main()
