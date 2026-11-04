"""Classify a blocked deltoid entry corridor before choosing any geometry cut.

Samples beyond the first collision are diagnostics, never reachable poses.
Any interference outside declared end regions or in R2 wire protection stops
a terminal-only relief interpretation. No source material is removed here.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,rod

DEST=A.DESIGN/'s1_deltoid_installation'


def main():
    base_file=DEST/'muscle_entry_study.json';base=read(base_file)
    file=ROOT/base['assembly_report'];r=read(file);assert sha(file)==base['assembly_sha256']
    solids={k:load_component(s) for k,s in r['sources'].items()};muscle=solids['tpu:deltoid']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==r['routes_sha256']
    own=[x for x in routes if x['head'] in ['FJ1468','FJ1513','FJ1467']]
    mf=A.I._manifold();zones=[];guards=[]
    for route in own:
        p=route['ports']['origin'];I=np.array(p['inner_face_source_mm']);n=A.D.unit(p['normal_source'])
        zones.append(rod(I-n*18.,I+n*.2,12.))
        guards.append(A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),2.))
    p=own[0]['ports']['distal'];I=np.array(p['inner_face_source_mm'])-np.array(p['offset_from_seat_mm']);n=A.D.unit(p['normal_source'])
    zones.append(rod(I-n*18.,I+n*.2,18.));zone=mf.Manifold.batch_boolean(zones,mf.OpType.Add)
    direction=A.D.unit([-1.,-.25,.25]);rows=[];first=None
    offsets=np.unique(np.r_[np.arange(0.,80.0001,1.),np.arange(0.,2.0001,.1)])
    fixed={k:solids[k] for k in base['fixed_keys']}
    for distance in offsets:
        hits=[]
        for key,solid in fixed.items():
            posed=muscle.translate(direction*distance);v=float(A.I.intersection_volume(posed,solid))
            if v<=.001:continue
            # Express overlap back in the unposed muscle's source frame.
            conflict=muscle^solid.translate(-direction*distance)
            hits.append(dict(target=key,volume_mm3=v,outside_declared_terminals_mm3=float((conflict-zone).volume()),
                protected_wire_R2_overlap_mm3=[dict(head=route['head'],volume_mm3=float(A.I.intersection_volume(conflict,guard)))
                    for route,guard in zip(own,guards)]))
        if hits and first is None:first=float(distance)
        rows.append(dict(offset_mm=float(distance),hits=hits,reachable_witness=False if first is not None else True))
        if abs(distance-round(distance))<1e-8:print('offset',distance,'targets',[h['target'] for h in hits],flush=True)
    result=dict(scope=__doc__,assembly_sha256=sha(file),muscle_entry_study_sha256=sha(base_file),routes_sha256=sha(route_file),
        script_sha256=sha(Path(__file__)),direction_source=direction.tolist(),first_collision_mm=first,
        sampled_corridor_length_mm=80.,fine_step_mm=.1,fine_interval_mm=[0.,2.],coarse_step_mm=1.,rows=rows,
        all_interference_inside_terminal_zones=not any(h['outside_declared_terminals_mm3']>.001 for row in rows for h in row['hits']),
        protected_wire_region_untouched_by_proposed_interference_cuts=not any(w['volume_mm3']>.001 for row in rows for h in row['hits'] for w in h['protected_wire_R2_overlap_mm3']),
        material_removed=False,source_geometry_unchanged=True,complete_assembly_verified=False,stage_pass=False)
    save(DEST/'entry_relief_analysis.json',result)
    print('terminal-only',result['all_interference_inside_terminal_zones'],'R2 safe',result['protected_wire_region_untouched_by_proposed_interference_cuts'],flush=True)


if __name__=='__main__':main()
