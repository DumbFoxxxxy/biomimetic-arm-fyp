"""Compare carrying four loose supports with the full deltoid during entry.

This is a candidate handling order, not a designed temporary fixture. Sources
are unchanged, no flexible bending is assumed, and bones stay present. A
clear rigid bundle still needs a real support-holding and staged fastening plan.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component

DEST=A.DESIGN/'s1_deltoid_installation'


def main():
    base_file=DEST/'muscle_entry_study.json';base=read(base_file)
    assert base['search_complete'] and not base['any_rigid_entry_clear_to_80mm']
    file=ROOT/base['assembly_report'];r=read(file);assert sha(file)==base['assembly_sha256']
    solids={k:load_component(s) for k,s in r['sources'].items()}
    moving=[base['moving_key']]+[k for k in base['required_early_keys'] if k.endswith('/front_support')]
    assert len(moving)==5
    fixed={k:solids[k] for k in base['fixed_keys'] if k not in moving}
    offsets=np.unique(np.r_[np.arange(0.,80.0001,1.),np.arange(0.,3.0001,.1)])
    rows=[]
    def record(complete):
        save(DEST/'support_bundle_entry_study.json',dict(scope=__doc__,assembly_report=file.relative_to(ROOT).as_posix(),
            assembly_sha256=sha(file),muscle_only_study_sha256=sha(base_file),script_sha256=sha(Path(__file__)),
            moving_keys=moving,fixed_keys=sorted(fixed),later_installation_keys=base['later_installation_keys'],
            final_inventory_preserved=True,search_complete=complete,directions_planned=len(base['rows']),directions_finished=len(rows),
            coarse_step_mm=1.,fine_interval_mm=[0.,3.],fine_step_mm=.1,corridor_length_mm=80.,intersection_threshold_mm3=.001,
            rows=rows,any_rigid_bundle_clear_to_80mm=any(x['clear_to_80mm'] for x in rows),
            temporary_support_holding_fixture_designed=False,whole_shoulder_installation_verified=False,TPU_deformation_verified=False,stage_pass=False))
    for source in base['rows']:
        direction=np.array(source['direction_source']);samples=[];blocked=None
        for offset in offsets:
            hits=[]
            for key in moving:
                posed=solids[key].translate(direction*offset)
                for target,solid in fixed.items():
                    v=float(A.I.intersection_volume(posed,solid))
                    if v>.001:hits.append(dict(moving=key,target=target,volume_mm3=v))
            samples.append(dict(offset_mm=float(offset),hits=hits))
            if hits:blocked=float(offset);break
        rows.append(dict(y_weight=source['y_weight'],z_weight=source['z_weight'],direction_source=direction.tolist(),
            samples=samples,first_collision_mm=blocked,clear_to_80mm=blocked is None))
        record(False)
        print(source['y_weight'],source['z_weight'],'first',blocked,'hits',samples[-1]['hits'][:3],flush=True)
    record(True)
    print('finished',len(rows),'clear',sum(x['clear_to_80mm'] for x in rows),flush=True)


if __name__=='__main__':main()
