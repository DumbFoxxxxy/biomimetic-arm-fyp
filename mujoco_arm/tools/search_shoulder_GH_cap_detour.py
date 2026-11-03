"""Finite two-leg route for cap 1, avoiding the humerus after axial release."""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from open_shoulder_GH_assembly_entries import DEST
from audit_shoulder_GH_assembly import scan
from check_shoulder_front_liner_keepers import dense_slide


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    joint=json.loads((A.G1/'assembly.json').read_text(encoding='utf8'));n=np.array(joint['neck_axis'])
    x=A.D.unit(np.cross(n,[0.,-1.,0.]));F=np.column_stack([x,np.cross(n,x),n])
    keys=['GH_cap_1','GH_top_washer_1']
    fixed_keys={k for k in S if k.startswith(('fixture:','bench:','clamp:','liner_keeper:'))
                or k in ['bone:scapula','bone:humerus','clavicle:fixed_support','GH_cap_0','GH_top_washer_0',
                          'GH_nut_0','GH_bottom_washer_0','GH_nut_1','GH_bottom_washer_1'] or k.startswith('clavicle:M4_')}
    fixed={k:S[k] for k in fixed_keys};attempts=[];selected=None
    directions=[A.D.unit(F@np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle)),slope]))
                for slope in [0.,.5,-.5] for angle in range(0,360,5)]
    for axial in [5.,10.,15.]:
        start={k:S[k].translate(n*axial) for k in keys}
        union=start[keys[0]]+start[keys[1]]
        for direction in directions:
            path=scan(union,direction,fixed,25.,.1)
            first=path['first_interference']
            attempts.append(dict(axial_release_mm=axial,direction_source=direction.tolist(),
                                 entire_sampled_path_clear=path['entire_sampled_path_clear'],
                                 first_hit_mm=None if first is None else first['distance_mm']))
            if path['entire_sampled_path_clear']:
                legs=[dense_slide({k:S[k] for k in keys},fixed,n,axial),
                      dense_slide(start,fixed,direction,25.)]
                if all(r['sampled_path_clear'] for r in legs):
                    selected=dict(axial_release_mm=axial,exit_direction_source=direction.tolist(),
                                  exit_translation_mm=25.,legs=legs,
                                  withdrawal_order=[0,1],placement_order=[1,0])
                    break
        if selected:break
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),moving_keys=keys,fixed_keys=sorted(fixed_keys),
                not_yet_installed_keys=sorted(set(S)-fixed_keys-set(keys)),attempts=attempts,selected=selected,
                stage_pass=False,actual_drive_range_deg=None)
    (DEST/'cap1_detour_search.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Attempts',len(attempts),'selected',None if selected is None else
          (selected['axial_release_mm'],selected['exit_direction_source']),flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
