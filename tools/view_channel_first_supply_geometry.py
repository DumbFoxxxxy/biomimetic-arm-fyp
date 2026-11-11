"""Rotate the neutral external-supply reference; zero-DOF scene is not a drive simulation."""
import argparse
import sys
import time
from build_s1_teres_origin_mount import ROOT,read
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_elbow_supply'
    record=read(folder/'neutral_scene.json');file=ROOT/record['model']['file'];assert provenance(file)==record['model']
    model=load_mujoco_model(file);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert model.nv==0 and model.ngeom==len(record['sources'])
    print('Neutral239-mesh reference only: green hollow supply rest shapes, original red elbow lattices, actual blue local sternum mount. No deformation or moving external-supply simulation. Camera can rotate/pan/zoom.',flush=True)
    if args.check:
        for source in record['sources'].values():assert provenance(ROOT/source['file'])['sha256']==source['sha256']
        print('Verified model and239 source fingerprints;0DOF',flush=True);return
    from mujoco import viewer
    with viewer.launch_passive(model,data) as handle:
        handle.cam.lookat[:]=[-.08,-.13,1.12];handle.cam.distance=.90;handle.cam.azimuth=55.;handle.cam.elevation=8.
        while handle.is_running():handle.sync();time.sleep(.01)


if __name__=='__main__':main()
