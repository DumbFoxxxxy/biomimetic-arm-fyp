"""Search zero-connected rigid placement paths for the actual flexible muscle.

A clear path is a nominal placement witness. A blocked rigid path does not
prove a flexible muscle impossible to install; deformation remains unverified.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component

DEST=A.DESIGN/'s1_teres_module_installation'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_teres_wire_stop/assembly_check.json';report=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['geometry_checks_pass']
    solids={key:load_component(row) for key,row in report['sources'].items()}
    # Parts to be installed after the bare muscle, not collision exclusions
    # in the final assembly. Record the full stage inventory explicitly.
    later={key for key in solids if key.startswith(('origin:FJ1507/','distal:FJ1507/','stop:FJ1507/','wire_active:FJ1507/'))}
    later|={'cap:FJ1507/origin','cap:FJ1507/distal'}
    moving='tpu:teres_major';fixed={key:solid for key,solid in solids.items() if key not in later|{moving}}
    muscle=solids[moving]
    base={key:float(A.I.intersection_volume(muscle,solid)) for key,solid in fixed.items()}
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    route=next(row for row in read(routes_file)['routes'] if row['head']=='FJ1507')
    a=np.array(route['ports']['origin']['normal_source']);b=np.array(route['ports']['distal']['normal_source'])
    bisector=A.D.unit(a+b);cross=A.D.unit(np.cross(a,b));difference=A.D.unit(a-b)
    directions=[]
    for c,d in [(0.,0.),(.25,0.),(-.25,0.),(.5,0.),(-.5,0.),(1.,0.),(-1.,0.),
                (0.,.15),(0.,-.15),(.5,.15),(-.5,.15),(.5,-.15),(-.5,-.15)]:
        direction=A.D.unit(bisector+cross*c+difference*d)
        if direction@a>0 and direction@b>0:directions.append((c,d,direction))
    rows=[]
    for c,d,direction in directions:
        samples=[];blocked=None
        for distance in np.arange(.5,40.0001,.5):
            posed=muscle.translate(direction*distance);hits=[]
            for key,solid in fixed.items():
                v=float(A.I.intersection_volume(posed,solid))
                if v>max(.001,base[key]+.001):hits.append(dict(key=key,volume_mm3=v,zero_pose_mm3=base[key]))
            samples.append(dict(offset_mm=float(distance),hits=hits))
            if hits:blocked=float(distance);break
        rows.append(dict(cross_weight=c,normal_difference_weight=d,direction_source=direction.tolist(),
            first_blocked_mm=blocked,last_sampled_no_increase_mm=(blocked-.5 if blocked else 40.),samples=samples,
            clear_to_40mm=blocked is None))
        print(c,d,'first',blocked,'hits',samples[-1]['hits'][:2],flush=True)
    save(DEST/'muscle_entry_study.json',dict(scope=__doc__,assembly_report=file.relative_to(ROOT).as_posix(),
        assembly_sha256=sha(file),stop_verification_sha256=sha(proof_file),script_sha256=sha(Path(__file__)),
        fixed_keys=sorted(fixed),later_installation_keys=sorted(later),moving_key=moving,
        no_part_deleted_from_final_assembly=True,step_mm=.5,threshold_increase_mm3=.001,
        inherited_zero_overlaps=[dict(key=k,volume_mm3=v) for k,v in base.items() if v>.001],
        rows=rows,TPU_deformation_verified=False,complete_shoulder_installation_verified=False,stage_pass=False))


if __name__=='__main__':main()
