"""Study continuous rigid placement of the complete actual deltoid lattice.

Muscle is tested rigidly as a sufficient placement witness only. Failure does
not prove TPU cannot bend into place. Preserve final inventory; every absent
part is explicitly later-installed and needs its own subsequent path proof.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component

DEST=A.DESIGN/'s1_deltoid_installation'
HEADS=['FJ1468','FJ1513','FJ1467']


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_deltoid_distal_mount/assembly_check.json';r=read(file)
    proof_file=file.parent/'whole_neutral_audit.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['neutral_pose_source_and_wire_pass']
    solids={k:load_component(s) for k,s in r['sources'].items()};moving='tpu:deltoid';muscle=solids[moving]
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==r['routes_sha256']
    normals=[A.D.unit(x['ports']['origin']['normal_source']) for x in routes if x['head'] in HEADS]
    normals.append(A.D.unit(next(x for x in routes if x['head']=='FJ1468')['ports']['distal']['normal_source']))
    later={k for k in solids if k.startswith(tuple(f'origin:{head}/' for head in HEADS)+('distal:FJ1468/',))}
    later|={k for k in solids if k.startswith(tuple(f'cap:{head}/origin' for head in HEADS))}
    later|={'cap:FJ1468/distal'}|{f'stop:{h}_stop_blank' for h in HEADS}
    required_early={f'origin:{h}/front_support' for h in HEADS}|{'distal:FJ1468/front_support','origin:FJ1467/nut_0','origin:FJ1467/nut_1'}
    assert required_early<=later
    later-=required_early
    fixed={k:s for k,s in solids.items() if k not in later|{moving}}
    rows=[]
    directions=[]
    for y in np.arange(-.75,.7501,.25):
        for z in np.arange(-.5,.7501,.25):
            direction=A.D.unit([-1.,y,z])
            if min(direction@n for n in normals)>.02:directions.append((float(y),float(z),direction))
    offsets=np.unique(np.r_[np.arange(0.,80.0001,1.),np.arange(0.,3.0001,.1)])
    def record(complete):
        save(DEST/'muscle_entry_study.json',dict(scope=__doc__,assembly_report=file.relative_to(ROOT).as_posix(),
            assembly_sha256=sha(file),neutral_audit_sha256=sha(proof_file),script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
            moving_key=moving,fixed_keys=sorted(fixed),required_early_keys=sorted(required_early),later_installation_keys=sorted(later),
            final_inventory_preserved=True,directions_planned=len(directions),directions_finished=len(rows),search_complete=complete,
            fine_interval_mm=[0.,3.],fine_step_mm=.1,coarse_step_mm=1.,corridor_length_mm=80.,intersection_threshold_mm3=.001,
            rows=rows,any_rigid_entry_clear_to_80mm=any(x['clear_to_80mm'] for x in rows),
            TPU_deformation_verified=False,whole_shoulder_installation_verified=False,stage_pass=False))
    for y,z,direction in directions:
        samples=[];blocked=None
        for offset in offsets:
            posed=muscle.translate(direction*offset);hits=[]
            for key,solid in fixed.items():
                v=float(A.I.intersection_volume(posed,solid))
                if v>.001:hits.append(dict(key=key,volume_mm3=v))
            samples.append(dict(offset_mm=float(offset),hits=hits))
            if hits:blocked=float(offset);break
        rows.append(dict(y_weight=y,z_weight=z,direction_source=direction.tolist(),samples=samples,
            first_collision_mm=blocked,clear_to_80mm=blocked is None))
        record(False)
        print(y,z,'first collision',blocked,'hits',samples[-1]['hits'][:3],flush=True)
    record(True)
    print('finished',len(rows),'clear',sum(x['clear_to_80mm'] for x in rows),flush=True)


if __name__=='__main__':main()
