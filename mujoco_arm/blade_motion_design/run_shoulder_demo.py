"""Play the actual eight-tendon shoulder; qpos is initialized once only."""
import argparse,json,time
import numpy as np
import mujoco
from shoulder_model import compile_model,inputs,EXPORT,v8,input_fingerprints

ACTIONS=['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--action',choices=ACTIONS,default='UP_ROT');parser.add_argument('--headless',action='store_true');parser.add_argument('--seconds',type=float,default=0.);args=parser.parse_args()
    report=json.loads((EXPORT/'shoulder_drive_report.json').read_text())
    if not report['pass'] or report['input_sha256']!=input_fingerprints():raise RuntimeError('Re-run shoulder drive audit for the current inputs')
    tables=np.load(EXPORT/'shoulder_control_tables.npz');model=compile_model();data=mujoco.MjData(model)
    cfg,fit,anchors=inputs();v8.set_joint_state(model,data,np.array(fit['poses'][0]['q_rad']));mujoco.mj_forward(model,data)
    choice={'action':args.action};current=args.action;cycle=28.;step=0;viewer=None
    if not args.headless:
        from mujoco import viewer as mjviewer
        def key_callback(key):
            if 49<=key<=54:choice['action']=ACTIONS[key-49]
        viewer=mjviewer.launch_passive(model,data,key_callback=key_callback)
        print('Keys 1-6 queue a motion for the next cycle: '+', '.join(ACTIONS))
    try:
        while (viewer is None or viewer.is_running()) and (not args.seconds or data.time<args.seconds):
            before=time.perf_counter();t=data.time%cycle
            if t<model.opt.timestep:current=choice['action']
            raw=t/12. if t<12. else (1. if t<14. else max(0.,1.-(t-14.)/12.))
            fraction=np.clip(raw,0,1);fraction=fraction*fraction*(3-2*fraction)
            controls=tables[f'{current}_ctrl'];position=fraction*(len(controls)-1);index=min(int(position),len(controls)-2);blend=position-index
            data.ctrl[:]=controls[index]*(1-blend)+controls[index+1]*blend
            mujoco.mj_step(model,data);step+=1
            if viewer is not None:
                if step%20==0:viewer.sync()
                time.sleep(max(0.,model.opt.timestep-(time.perf_counter()-before)))
        print({'finite':bool(np.isfinite(data.qpos).all()),'simulated_seconds':float(data.time),'steps':step})
    finally:
        if viewer is not None:viewer.close()
if __name__=='__main__':main()
