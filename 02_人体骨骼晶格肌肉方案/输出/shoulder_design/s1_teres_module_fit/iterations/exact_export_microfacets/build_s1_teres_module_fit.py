"""Local terminal revision for muscle entry and positive distal screw capture.

Working lattice and bone remain unchanged. Sampled insertion relief and a TPU
flange are geometric interfaces, not strength, preload or driven-motion proof.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,shoulder_pads,sha,read,save,export,load_component

DEST=A.DESIGN/'s1_teres_module_fit'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_teres_wire_stop/assembly_check.json';parent=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['geometry_checks_pass']
    entry_file=A.DESIGN/'s1_teres_module_installation/muscle_entry_study.json';entry=read(entry_file)
    assert entry['assembly_sha256']==sha(file)
    best=min(entry['rows'],key=lambda row:sum(hit['volume_mm3'] for hit in row['samples'][0]['hits']))
    direction=np.array(best['direction_source'])
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(route_file)['routes'];route=next(row for row in routes if row['head']=='FJ1507')
    source=parent['sources']['tpu:teres_major'];muscle=load_component(source)
    bone_row=parent['sources']['bone:scapula'];bone=load_component(bone_row)
    origin=np.array(route['ports']['origin']['inner_face_source_mm']);n0=A.D.unit(route['ports']['origin']['normal_source'])
    terminal_zone=rod(origin-n0*20.,origin+n0*.2,12.)
    distal=np.array(route['ports']['distal']['inner_face_source_mm']);n1=A.D.unit(route['ports']['distal']['normal_source'])
    F=shoulder_pads.pad_basis(n1)
    flange=rod(distal-n1*3.2,distal-n1*1.,10.5)
    for sign in [-1,1]:
        centre=distal+F[:,0]*(sign*7)
        flange=flange-rod(centre-n1*4.,centre+n1*1.,1.7)
    for wire in routes:
        flange=flange-A.D.outer_envelope(np.array(wire['full_wire_path_source_mm']),.86)
    revised=muscle+flange
    trace=[];outside=[]
    distances=np.unique(np.r_[np.arange(0.,40.0001,.5),np.arange(0.,2.0001,.1)])
    # Subtract only actual interfering bone stock in the declared proximal
    # terminal. Never trim working cells or another muscle to pass a path.
    for distance in distances:
        conflicting=muscle^bone.translate(-direction*distance)
        volume=float(conflicting.volume())
        if volume>1e-9:
            beyond=float((conflicting-terminal_zone).volume())
            trace.append(dict(offset_mm=float(distance),original_overlap_mm3=volume,outside_zone_mm3=beyond))
            if beyond>.001:outside.append(trace[-1])
            revised=revised-(conflicting^terminal_zone)
    actual,meta=export(revised,DEST/'teres_screw_flange_and_lead_in.obj',np.array(source['translation_source_mm']))
    sources=dict(parent['sources']);sources['tpu:teres_major']=meta
    solids={key:(actual if key=='tpu:teres_major' else load_component(row)) for key,row in sources.items()}
    hits=[]
    for key,solid in solids.items():
        if key=='tpu:teres_major':continue
        volume=float(A.I.intersection_volume(actual,solid))
        before=float(A.I.intersection_volume(muscle,solid))
        if volume>.001 or before>.001:hits.append(dict(key=key,volume_mm3=volume,parent_mm3=before,
            new_or_increased=volume>max(.001,before+.001)))
    removed=muscle-actual;added=actual-muscle
    distal_zone=rod(distal-n1*4.,distal+n1*.2,12.)
    wall_guard=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),2.)
    guard_loss=float((removed^wall_guard).volume())
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),
        parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),
        entry_study_sha256=sha(entry_file),script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        replaced_keys=['tpu:teres_major'],new_keys=[],all_other_206_sources_unchanged=True,
        entry_direction_source=direction.tolist(),entry_trace=trace,unresolved_outside_terminal=outside,
        dimensions=dict(proximal_relief_radius_mm=12.,proximal_relief_normal_mm=[-20.,.2],
            relief_coarse_step_mm=.5,relief_fine_interval_mm=[0.,2.],relief_fine_step_mm=.1,
            distal_flange_radius_mm=10.5,distal_flange_normal_mm=[-3.2,-1.],screw_pitch_mm=14.,flange_holes_mm=3.4),
        muscle_removed_mm3=float(removed.volume()),muscle_added_mm3=float(added.volume()),
        removed_outside_proximal_terminal_mm3=float((removed-terminal_zone).volume()),
        added_outside_distal_terminal_mm3=float((added-distal_zone).volume()),
        removed_within_wire_outer_2mm_guard_mm3=guard_loss,changed_part_pairs=206,material_intersections=hits,
        no_new_or_increased_overlap=not any(row['new_or_increased'] for row in hits),
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners=parent['intended_body_owners'],adopted_as_current=False,stage_pass=False,
        actual_drive_range_deg=None,TPU_strength_verified=False,complete_installation_verified=False)
    save(DEST/'assembly_check.json',report)
    print({k:report[k] for k in ['muscle_removed_mm3','muscle_added_mm3','removed_outside_proximal_terminal_mm3',
        'removed_within_wire_outer_2mm_guard_mm3','no_new_or_increased_overlap']},flush=True)
    print('Outside-zone trace',outside,'changed pair hits',hits,flush=True)
    if outside or guard_loss>.001 or not report['no_new_or_increased_overlap']:
        raise ValueError('Terminal-only revision does not resolve this candidate path; saved evidence retained')


if __name__=='__main__':main()
