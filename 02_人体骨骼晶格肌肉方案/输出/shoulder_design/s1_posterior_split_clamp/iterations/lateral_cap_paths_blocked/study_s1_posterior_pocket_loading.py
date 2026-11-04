"""Evaluate loading nuts after caps and preserving the B nut-bearing floor.

This proposes metal entry relief only. Actual modified registration and final
geometry need generation and independent replay before any local pass claim.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,rod


def main():
    folder=A.DESIGN/'s1_posterior_split_clamp';file=folder/'assembly_check.json';r=read(file)
    solids={key:load_component(row) for key,row in r['sources'].items()}
    I=np.array(r['local_origin_source_mm']);F=np.array(r['local_basis_columns']);side,t,n=F.T
    nuts=[]
    for index in range(2):
        key=f'origin:FJ1513/nut_{index}';hits=[]
        for offset in np.arange(0.,3.0001,.1):
            posed=solids[key].translate(n*offset)
            for target,obstacle in solids.items():
                if target==key:continue
                v=float(A.I.intersection_volume(posed,obstacle))
                if v>.001:hits.append(dict(offset_mm=float(offset),target=target,volume_mm3=v))
        nuts.append(dict(index=index,hits=hits))
    cap=solids['cap:FJ1513/origin_B'];revised=cap
    # Clip only the obstacle search volume; every cap stock change is measured.
    box=A.I._manifold().Manifold.cube([60.,40.,8.]).translate([-30.,-20.,-1.]).transform(np.c_[F,I])
    obstacles=[solids[key]^box for key in ['tpu:deltoid','tpu:infraspinatus']]
    for distance in np.unique(np.r_[np.arange(0.,12.0001,.5),np.arange(0.,3.0001,.1)]):
        for obstacle in obstacles:revised=revised-obstacle.translate(side*distance)
    parent=read(ROOT/r['parent_report']);index=int(next(k for k,v in r['nut_to_half'].items() if v=='B'))
    p=np.array(parent['selected_bone_crossings'][index]['axis_source_mm'])
    witness=rod(p+n*.3,p+n*1.7,2.6)-rod(p+n*.2,p+n*1.8,1.75)
    result=dict(source_sha256=sha(file),script_sha256=sha(Path(__file__)),nut_axial_entry_mm=3.,nut_step_mm=.1,
        nuts=nuts,B_opening_direction_source=(-side).tolist(),B_cap_entry_relief_removed_mm3=float((cap-revised).volume()),
        B_nut_floor_witness_missing_mm3=float((witness-revised).volume()),stage_pass=False)
    save(folder/'pocket_loading_study.json',result)
    print(dict(nut_hits=[len(row['hits']) for row in nuts],
        cap_removed=result['B_cap_entry_relief_removed_mm3'],floor_missing=result['B_nut_floor_witness_missing_mm3']),flush=True)


if __name__=='__main__':main()
