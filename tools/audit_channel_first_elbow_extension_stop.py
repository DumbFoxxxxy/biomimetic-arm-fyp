"""Solid stop geometry sweep; prescribed poses are clearance checks, not actuation."""
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    assembly_file=base/'elbow_extension_stop/assembly.json';assembly=read(assembly_file)
    check_file=base/'elbow_passive_model/model_check.json';check=read(check_file)
    owners={r['source_key']:r['owner'] for r in check['geom_owners']}
    owners.update({k:'humerus' for k in assembly['sources'] if k.startswith('stop:')})
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    E=np.array(assembly['origin_source_mm']);F=np.array(assembly['frame_columns'])
    moving={k for k,v in owners.items() if v in ('forearm','FJ1486_distal_pivot','FJ1477_distal_pivot')}
    pairs=[(a,b) for a,b in combinations(solids,2) if (a in moving)!=(b in moving)]
    def hits(angle):
        rotation=A.I.rotation_about(E,F[:,0],angle)
        placed={k:s.transform(rotation) if k in moving else s for k,s in solids.items()}
        result=[]
        for a,b in pairs:
            volume=float(A.I.intersection_volume(placed[a],placed[b]))
            if volume>.001:result.append(dict(left=a,right=b,overlap_mm3=volume))
        return result
    zero=hits(0.);assert not zero,zero
    rows=[];last=0.;first=None
    for angle in range(1,31):
        collision=hits(angle);rows.append(dict(angle_deg=angle,hits=collision))
        if collision:
            lo,hi=last,float(angle)
            while hi-lo>.10000001:
                mid=(lo+hi)/2
                if hits(mid):hi=mid
                else:lo=mid
            first=dict(clear_deg=lo,blocked_deg=hi,hits=hits(hi));break
        last=float(angle)
    # Check the return poses explicitly, with cap swivels held at their neutral angles.
    returning=[dict(angle_deg=a,hits=hits(a)) for a in np.arange(last,-.1,-1.)]
    negative=[dict(angle_deg=a,hits=hits(a)) for a in (-.1,-.2,-.3,-.4)]
    save(assembly_file.parent/'motion_screen.json',dict(assembly=provenance(assembly_file),
        owner_reference=provenance(check_file),generator=provenance(ROOT/'tools/audit_channel_first_elbow_extension_stop.py'),
        rigid_parts=len(solids),varying_pairs=len(pairs),zero_hits=zero,outgoing=rows,returning=returning,
        first_positive_interference=first,negative_stop_probes=negative,sampling_deg=1.,boundary_resolution_deg=.1,
        overlap_threshold_mm3=.001,
        scope='Rigid geometric outgoing/return elbow-only screening through at most 30deg. Both cap hinges remain neutral relative to forearm; no force-drive or full K-range acceptance. All changing pairs included, no whole pair exclusions. Zero-gap plane contact is allowed only as zero-volume contact. Same-body assembly fit checked by generator. No contact force, tool access or insertion sweep validation.'))
    print('Last clear integer pose',last,'first positive interference',first,'negative .1deg',negative[0],flush=True)


if __name__=='__main__':main()
