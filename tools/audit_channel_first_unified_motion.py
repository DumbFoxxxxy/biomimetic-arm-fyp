"""Sample actual solids under ideal joint poses; not force-driven acceptance."""
from itertools import combinations
import argparse
import sys

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    import mujoco
    sys.path.insert(0, str(ROOT))
    from project_runtime import load_mujoco_model
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_skeleton'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder')
    args=parser.parse_args()
    if args.folder:folder=ROOT/args.folder
    check_file = folder/'geometry_mass_body_chain_check.json'
    check = read(check_file)
    model_file = ROOT/check['generated_model']['file']
    assert provenance(model_file) == check['generated_model']
    reference_file = ROOT/check['reference']['file']
    assert provenance(reference_file) == check['reference']
    reference = read(reference_file)
    sources = {**reference['sources'], **reference['joint_context_sources']}
    solids = {key: load_component(value) for key, value in sources.items()}
    owners = {r['source_key']: r['owner'] for r in check['geom_owners']}
    pairs = [(a,b) for a,b in combinations(solids,2) if owners[a] != owners[b]]
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model,data)
    origins = {name: np.array(value)*1000. for name,value in check['neutral_body_origins_m'].items()}
    origins['thorax'] = np.zeros(3)
    body_ids = {name: mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,name) for name in origins}
    def hits(delta):
        mujoco.mj_resetData(model,data)
        mujoco.mj_integratePos(model,data.qpos,delta,1.)
        mujoco.mj_forward(model,data)
        transforms = {}
        for name,i in body_ids.items():
            rotation = data.xmat[i].reshape(3,3)
            transforms[name] = np.c_[rotation, data.xpos[i]*1000.-rotation@origins[name]]
        placed = {key: solid.transform(transforms[owners[key]]) for key,solid in solids.items()}
        result = []
        for a,b in pairs:
            volume = float(A.I.intersection_volume(placed[a],placed[b]))
            if volume > .001:
                result.append(dict(left=a,right=b,overlap_mm3=volume))
        return result
    zero = hits(np.zeros(model.nv))
    rows = []
    if not zero:
        for joint,axis,sign,limit in [(joint,a,s,10) for joint in ('SC','AC','GH') for a in range(3) for s in (-1,1)]+[('elbow',0,1,30)]:
            j = mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,joint)
            dof = model.jnt_dofadr[j]+axis
            last, first = 0., None
            for angle in range(1,limit+1):
                delta = np.zeros(model.nv);delta[dof]=np.radians(sign*angle)
                collision = hits(delta)
                if collision:
                    first = dict(angle_deg=sign*angle,hits=collision)
                    break
                last = sign*angle
            rows.append(dict(joint=joint,local_axis=axis,direction=sign,
                             last_clear_sample_deg=last,first_interference=first))
            print(rows[-1],flush=True)
    save(folder/'motion_screen.json', dict(model=provenance(model_file),
        body_chain_check=provenance(check_file), reference=provenance(reference_file),
        tested_cross_body_pairs=len(pairs), baseline_interferences=zero, small_angle_probes=rows,
        sampling_deg=1., threshold_mm3=.001,
        status='baseline_interferences_require_resolution' if zero else 'sampled_small_angle_only',
        scope='All cross-rigid-body solid pairs included, including adjacent bones and joint hardware. Same-body pairs do not change under these poses; internal assembly fit remains a separate task. No licensed contact regions assumed. No K-range or force-drive acceptance. Boundaries not yet refined; paths are outgoing probes only, not final outward/return validation.'))
    print('Baseline hits',zero,'cross-body pairs',len(pairs),flush=True)


if __name__ == '__main__':
    main()
