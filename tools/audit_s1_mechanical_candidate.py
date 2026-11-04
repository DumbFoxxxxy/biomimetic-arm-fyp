"""Full zero-pose source audit for a mechanical working shoulder candidate.

Every distinct source pair is tested; no adjacent-bone/contact pair is omitted.
Own wire reference geometry alone is omitted from its identical route probe.
This is a sampled neutral-pose audit, not motion, assembly, FE or physical proof.
"""
from pathlib import Path
import argparse
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',default=str((A.DESIGN/'s1_deltoid_distal_mount/assembly_check.json').relative_to(ROOT)))
    args=parser.parse_args();file=ROOT/args.report;r=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['static_fit_contact_bearing_pass'] and proof['local_approaches_geometry_pass']
    solids={k:load_component(s) for k,s in r['sources'].items()}
    for source in r['sources'].values():assert sha(ROOT/source['file'])==source['sha256']
    hits=[];count=0
    for a,b in combinations(solids,2):
        v=float(A.I.intersection_volume(solids[a],solids[b]));count+=1
        if v>.001:hits.append(dict(a=a,b=b,volume_mm3=v,includes_soft_geometry=a.startswith('tpu:') or b.startswith('tpu:')))
    assert count==len(solids)*(len(solids)-1)//2
    print('all pairs',count,'hits',hits,flush=True)
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==r['routes_sha256'];wires=[]
    for route in routes:
        head=route['head'];corridor=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3)
        wh=[];own_references=[]
        for key,solid in solids.items():
            if key==f'wire_active:{head}/reference':
                own_references.append(key);continue
            v=float(A.I.intersection_volume(corridor,solid))
            if v>.001:wh.append(dict(key=key,volume_mm3=v))
        wires.append(dict(head=head,distinct_source_probes=len(solids)-len(own_references),identical_wire_self_references=own_references,hits=wh))
        print('wire',head,'hits',wh,flush=True)
    missing=set(solids)-set(r['intended_body_owners']);stale=set(r['intended_body_owners'])-set(solids)
    result=dict(scope=__doc__,assembly_report=file.relative_to(ROOT).as_posix(),assembly_sha256=sha(file),
        local_verification_sha256=sha(proof_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        sources_checked=len(solids),all_source_hashes_match=True,neutral_pose_pairs=count,excluded_source_pairs=[],
        intersection_volume_threshold_mm3=.001,neutral_pose_intersections=hits,six_active_wire_corridors=wires,
        missing_body_owners=sorted(missing),stale_body_owners=sorted(stale),
        neutral_pose_source_and_wire_pass=not hits and not any(x['hits'] for x in wires) and not missing and not stale,
        original_stop_blanks_are_unverified_references=proof['unverified_late_stop_blank_keys'],
        complete_assembly_verified=False,driven_motion_verified=False,physical_hardware_strength_verified=False,
        S1_frozen=False,stage_pass=False,actual_drive_range_deg=None)
    save(file.parent/'whole_neutral_audit.json',result)
    print('neutral pass',result['neutral_pose_source_and_wire_pass'],flush=True)


if __name__=='__main__':main()
