"""Unilateral endpoint-force allocation on the actual new SC/AC mass model."""
import sys
import numpy as np
from scipy.optimize import lsq_linear

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance
from simulate_channel_first_coordinated_spools import endpoint_moment_matrix


def main():
    import mujoco
    sys.path.insert(0,str(ROOT))
    from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/unified_SC_AC'
    endpoint_file=folder/'arm_girdle_equivalent_drive_geometry_mass_check.json'
    endpoint=read(endpoint_file);model_file=ROOT/endpoint['generated_model']['file']
    assert provenance(model_file)==endpoint['generated_model']
    path_file=folder/'SC_AC_path_screen.json';paths=read(path_file)
    chain=read(ROOT/endpoint['chain_check']['file'])
    assert paths['chain_check']==endpoint['chain_check']
    assert provenance(ROOT/chain['reference']['file'])==chain['reference']
    model=load_mujoco_model(model_file);data=mujoco.MjData(model)
    blocks={k:v['velocity_dof_indices'] for k,v in endpoint['neutral_joint_blocks'].items()}
    rows=[]
    # Study all six geometrically selected branches, keeping boundary warnings.
    for axis in range(3):
        for sign in (-1,1):
            candidates=[r for r in paths['paths'] if r['axis']==axis and r['direction']==sign]
            chosen=max(candidates,key=lambda r:(abs(r['last_contiguous_clear_SC_deg']),-abs(r['endpoint']['AC_angle_deg'])))
            end=chosen['last_contiguous_clear_SC_deg'];ratio=chosen['AC_to_SC_ratio']
            steps=max(1,int(np.ceil(abs(end)*max(1.,abs(ratio)))))
            samples=[]
            for sc in np.linspace(0.,end,steps+1):
                mujoco.mj_resetData(model,data);delta=np.zeros(model.nv)
                delta[blocks['SC'][axis]]=np.radians(sc)
                delta[blocks['AC'][axis]]=np.radians(sc*ratio)
                mujoco.mj_integratePos(model,data.qpos,delta,1.);mujoco.mj_forward(model,data)
                matrix=endpoint_moment_matrix(model,data,endpoint['channels'])
                target=data.qfrc_bias-data.qfrc_passive
                fit=lsq_linear(matrix,target,bounds=(0.,60.),tol=1e-11,max_iter=1000)
                unbounded=lsq_linear(matrix,target,bounds=(0.,np.inf),tol=1e-11,max_iter=1000)
                samples.append(dict(SC_angle_deg=float(sc),AC_angle_deg=float(sc*ratio),
                    tensions_N=fit.x.tolist(),residual_Nm=fit.fun.tolist(),
                    residual_norm_Nm=float(np.linalg.norm(fit.fun)),peak_tension_N=float(fit.x.max()),
                    residual_by_joint_Nm={k:float(np.linalg.norm(fit.fun[v])) for k,v in blocks.items()},
                    no_upper_tension_bound_residual_Nm=float(np.linalg.norm(unbounded.fun)),
                    span_lengths_mm=(data.ten_length*1000).tolist()))
            rows.append(dict(axis=axis,direction=sign,ratio=ratio,
                endpoint_gap_mm=chosen['endpoint']['minimum_gap_mm_capped0_5'],
                actual_humerus_direction_change_deg=chosen['endpoint']['humerus_direction_change_deg'],
                warning='Near-contact endpoint is not a working pose' if chosen['endpoint']['minimum_gap_mm_capped0_5']<.1 else 'Geometry study limit only',
                samples=samples))
            print(axis,sign,'max residual',max(r['residual_norm_Nm'] for r in samples),
                  'unbounded at end',samples[-1]['no_upper_tension_bound_residual_Nm'],flush=True)
    save(folder/'endpoint_force_balance_paths.json',dict(
        endpoint_model=provenance(endpoint_file),geometric_paths=provenance(path_file),
        generator=provenance(ROOT/'tools/study_channel_first_unified_force_balance.py'),
        cases=rows,channels=[c['id'] for c in endpoint['channels']],tension_bounds_N=[0.,60.],
        scope='Static inverse force allocation only on new uniform-density assumed-mass skeleton. Two-point tendon surrogates, not actual cable routing. No passive scapulothoracic support, guide reaction, steel/conduit friction, hardware force or driven trajectory acceptance. Unlimited-tension residual tests positive-force direction, not a recommendation to increase pull.'))


if __name__=='__main__':main()
