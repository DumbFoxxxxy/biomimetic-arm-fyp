"""Screen oblique approaches of split caps with manually held front nuts.

Registration counterpart is omitted only during screening; final pins must
be aligned and independently replayed with both cap halves present.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component


def main():
    folder=A.DESIGN/'s1_posterior_split_clamp';file=folder/'iterations/lateral_cap_paths_blocked/assembly_check.json'
    r=read(file);proof=read(file.parent/'verification_checks.json')
    assert proof['assembly_sha256']==sha(file)
    solids={}
    for key,row in r['sources'].items():
        frozen=dict(row)
        if key in r['new_keys']+r['replaced_keys']:
            frozen['file']=(file.parent/Path(row['file']).name).relative_to(ROOT).as_posix()
        solids[key]=load_component(frozen)
    F=np.array(r['local_basis_columns']);side,tangent,n=F.T;groups=[]
    for suffix,sign in [('A',-1),('B',1)]:
        old_stage=next(row for row in proof['local_stages'] if row['name']=='cap_and_held_nut_'+suffix)
        fixed=set(old_stage['fixed_keys'])-{'cap:FJ1513/origin_A','cap:FJ1513/origin_B'}
        nut_by_half={half:'origin:FJ1513/nut_'+index for index,half in r['nut_to_half'].items()}
        if suffix=='B':fixed.discard(nut_by_half['A'])  # B first, A not yet installed.
        else:fixed.add(nut_by_half['B'])  # A second; B registration geometry must be regenerated.
        moving=old_stage['moving'];rows=[]
        cases=[(sw,tw) for sw in ([-1.] if suffix=='A' else [-1.,1.]) for tw in [-2.,-1.,-.5,0.,.5,1.,2.]]
        for side_weight,t_weight in cases:
            for n_weight in [0.,.25,.5,1.]:
                direction=A.D.unit(side*side_weight+tangent*t_weight+n*n_weight)
                nut_hits=[];cap_hits=[]
                for distance in [0.,.1,.2,.5,1.,2.,3.,4.,6.,8.,10.,12.]:
                    for key in moving:
                        posed=solids[key].translate(direction*distance)
                        for target in fixed:
                            v=float(A.I.intersection_volume(posed,solids[target]))
                            if v>.001:
                                hit=dict(offset_mm=distance,moving=key,target=target,volume_mm3=v)
                                (cap_hits if key.startswith('cap:') else nut_hits).append(hit)
                rows.append(dict(side_weight=side_weight,tangent_weight=t_weight,normal_weight=n_weight,direction_source=direction.tolist(),
                    nut_hits=nut_hits,cap_hits=cap_hits,coarse_clear=not cap_hits and not nut_hits))
                print(suffix,side_weight,t_weight,n_weight,'nut',len(nut_hits),'cap',len(cap_hits),flush=True)
        groups.append(dict(half=suffix,fixed_keys=sorted(fixed),moving=moving,rows=rows))
    save(folder/'cap_direction_study.json',dict(scope=__doc__,source_report=file.relative_to(ROOT).as_posix(),
        source_sha256=sha(file),source_verification_sha256=sha(file.parent/'verification_checks.json'),
        script_sha256=sha(Path(__file__)),archived_modified_meshes_used=True,groups=groups,
        assembly_order_screened=['B','A'],registration_alignment_pending=True,stage_pass=False))


if __name__=='__main__':main()
