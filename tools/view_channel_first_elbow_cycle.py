"""Interactive replay of calculated elastic-channel cycles, with audit status."""
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
    parser.add_argument('--trial',choices=['k1e6','k1e7','k1e8','buffer_k1e6'],default='buffer_k1e6')
    parser.add_argument('--speed',type=float,default=1.)
    parser.add_argument('--folder',help='Alternative model/result directory, relative to project root')
    parser.add_argument('--result-name',help='Alternative response basename, e.g. local_cycle')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args();assert args.speed>0
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_extension_stop'
    if args.folder:folder=ROOT/args.folder
    result_name=args.result_name or f'cycle_{args.trial}'
    assert all(c.isalnum() or c in '_-' for c in result_name)
    response=folder/(result_name+'.json')
    if not response.exists() and (folder/(result_name+'_numerical_rejection.json')).exists():
        print('Numerically invalid response retained in a verified archive; not replayed as mechanical motion.',flush=True);return
    r=read(response)
    model_file=ROOT/r['model']['file'];assert provenance(model_file)==r['model']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model);neutral=data.qpos.copy()
    runtime=r['runtime_contact_parameters'];model.pair_solref[:]=np.array(runtime['pair_solref'])
    if 'pair_margin_m' in runtime:model.pair_margin[:]=np.array(runtime['pair_margin_m'])
    trace=r['trace'];audit_file=folder/(result_name+'_solid_screen.json')
    if audit_file.exists():
        audit=read(audit_file);assert audit['response']==provenance(response)
        if audit.get('status')=='rejected_numerical_response':
            print('Numerically invalid response; not replayed as a mechanical motion:',audit['rejection_reasons'],flush=True);return
        first=audit['first_rigid_interference']
        for value in audit['first_interferences'].values():
            if value and (first is None or value['time_s']<first['time_s']):first=value
        if first:
            trace=[s for s in trace if s['time_s']<first['time_s']]
            print('Replay stops before first audited interference at',first['time_s'],'s; later poses rejected.',flush=True)
        else:print('Saved rigid and steel geometry audit clear; flexible wall, numerical contact and full routes still require review.',flush=True)
        if any(v is not None for v in audit.get('external_wall_first_intersections',{}).values()):
            print('Unsolved external-wall contact exists. Replay is a free-space diagnostic beyond that contact, not verified physical motion.',flush=True)
    else:print('Solid audit pending: this is an exploratory numerical replay.',flush=True)
    assert trace
    times=np.array([s['time_s'] for s in trace])
    def frame(index):
        state=trace[index];data.qpos[:]=neutral
        mujoco.mj_integratePos(model,data.qpos,np.array(state['generalized_displacement']),1.)
        data.qvel[:]=0.;data.ctrl[:]=state['tensions_N'];mujoco.mj_forward(model,data)
        assert np.isfinite(data.qpos).all()
    print('Whole thorax/arm, real stop and cap solids, two moving elbow channel beams. Other25 muscles are lines; full TPU muscle shapes and supply sheaths absent. Recorded motion, not a live new drive simulation.',flush=True)
    if (folder/'assembly.json').exists() and 'interfaces' in read(folder/'assembly.json'):
        print('Includes real supply-end hardware. Green short sheath meshes are neutral reference displays, not simulated chest supply guides.',flush=True)
    if args.check:
        for index in range(len(trace)):frame(index)
        print('Reconstructed',len(trace),'poses;',model.nv,'DOFs;',model.ngeom,'geoms',flush=True);return
    from mujoco import viewer
    with viewer.launch_passive(model,data) as v:
        v.cam.lookat[:]=[-.08,-.08,1.18];v.cam.distance=.85;v.cam.azimuth=125;v.cam.elevation=8
        start=time.perf_counter()
        while v.is_running():
            t=((time.perf_counter()-start)*args.speed)%float(times[-1])
            i=max(0,min(len(trace)-1,int(np.searchsorted(times,t,side='right')-1)))
            with v.lock():frame(i)
            v.sync();time.sleep(.01)


if __name__=='__main__':main()
