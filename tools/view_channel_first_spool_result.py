"""Replay computed force-driven states; this viewer does not validate new motion."""
import argparse
import sys
import time
import numpy as np
from build_s1_teres_origin_mount import ROOT,read
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial',choices=['hold','flexor_cycle','flexor_release'],default='flexor_cycle')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    file=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC/spool_drive_dynamic_screen.json'
    report=read(file);model_file=ROOT/report['model']['file']
    assert provenance(model_file)==report['model']
    audit=read(file.parent/'spool_dynamic_solid_audit.json')
    assert audit['source']==provenance(file)
    trial=next(t for t in report['trials'] if t['trial']==args.trial)
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    neutral=data.qpos.copy();trace=trial['trace'];times=np.array([r['time_s'] for r in trace])
    elbow=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'elbow')
    def frame(index):
        state=trace[index];data.qpos[:]=neutral
        mujoco.mj_integratePos(model,data.qpos,np.array(state['joint_displacement_rad']),1.)
        data.qvel[:]=0.;data.ctrl[:]=state['cable_tension_N'];mujoco.mj_forward(model,data)
        assert abs(np.degrees(data.qpos[model.jnt_qposadr[elbow]])-state['elbow_deg'])<1e-7
    print('Computed dynamics replay: endpoint cable proxies and pad envelopes; no manufactured TPU muscles or routed supply. Camera is interactive; motion is replayed, not recomputed here.',flush=True)
    if args.check:
        for index in range(len(trace)):frame(index)
        print('Verified',len(trace),'reconstructed states;',model.ngeom,'visual geoms',flush=True)
        return
    from mujoco import viewer
    with viewer.launch_passive(model,data) as v:
        v.cam.lookat[:]=[-.08,-.08,1.18];v.cam.distance=.85;v.cam.azimuth=125;v.cam.elevation=8
        start=time.perf_counter()
        while v.is_running():
            t=(time.perf_counter()-start)%report['duration_s']
            index=max(0,min(len(trace)-1,int(np.searchsorted(times,t,side='right')-1)))
            with v.lock():frame(index)
            v.sync();time.sleep(.01)


if __name__=='__main__':main()
