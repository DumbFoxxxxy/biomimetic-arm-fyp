"""Local distal TPU lead-in for staged teres installation with collar screws later.

Bone, metal carrier, working lattice and original wire coordinates are retained.
This prepares nominal geometry, not actual clamp grip or physical installation.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,sha,read,save,export,load_component,shoulder_pads

DEST=A.DESIGN/'s1_teres_assembly_order'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_teres_split_cap/assembly_check.json';parent=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['split_cap_geometry_pass']
    entry_file=A.DESIGN/'s1_teres_module_fit/assembly_check.json';entry=read(entry_file)
    assert parent['sources']['tpu:teres_major']==entry['sources']['tpu:teres_major']
    direction=np.array(entry['entry_direction_source'])
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    route=next(row for row in read(route_file)['routes'] if row['head']=='FJ1507')
    I=np.array(route['ports']['distal']['inner_face_source_mm']);n=A.D.unit(route['ports']['distal']['normal_source'])
    zone=rod(I-n*4.,I+n*.2,12.)
    muscle=load_component(parent['sources']['tpu:teres_major'])
    protected=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),2.)
    # Preserve the complete existing central bearing boss as well as the
    # wire wall. Its intentional support contact is not a request for an
    # artificial uniform 0.03mm gap across the bearing surface.
    protected_for_cut=protected+rod(I-n*4.1,I+n*.3,5.2)
    carrier=load_component(parent['sources']['collar:jaw_lower']);guard=carrier
    for axis in range(3):
        for sign in [-1,1]:guard=guard+carrier.translate(np.eye(3)[axis]*sign*.03)
    revised=muscle;revised_carrier=carrier;trace=[];outside=[]
    terminal_wall=muscle^protected_for_cut^zone
    metal_zone=rod(I-n*7.,I+n*.5,6.)
    F=shoulder_pads.pad_basis(n)
    thread_stock=None
    for sign in [-1,1]:
        centre=I+F[:,0]*(sign*7)
        witness=rod(centre-n*6.2,centre-n*3.2,2.8)-rod(centre-n*6.3,centre-n*3.1,1.61)
        thread_stock=witness if thread_stock is None else thread_stock+witness
    for offset in np.unique(np.r_[np.arange(0.,40.0001,.5),np.arange(0.,10.0001,.1)]):
        # Test actual carrier material. Its flat support face is an intended
        # contact, so a uniform inflation there would incorrectly shave the
        # bearing boss even where no insertion interference exists.
        shifted=carrier.translate(-direction*offset);conflict=muscle^shifted
        volume=float(conflict.volume())
        if volume<=1e-9:continue
        actual_outside=float(((muscle^carrier.translate(-direction*offset))-zone).volume())
        trace.append(dict(offset_mm=float(offset),guard_overlap_mm3=volume,
            actual_carrier_outside_terminal_mm3=actual_outside))
        if actual_outside>.001:outside.append(trace[-1])
        revised=revised-(conflict^zone)
        # Keep the TPU guide wall. Remove only the carrier material occupying
        # its swept guide passage, within a named metal region and outside
        # the two nominal M3 thread stock witnesses.
        moved_wall=terminal_wall.translate(direction*offset)
        revised_carrier=revised_carrier-((moved_wall^metal_zone)-thread_stock)
    prepared=revised.as_original().simplify(1e-6)
    actual,meta=export(prepared,DEST/'teres_bridge_lead_in.obj',np.array(parent['sources']['tpu:teres_major']['translation_source_mm']))
    sources=dict(parent['sources']);sources['tpu:teres_major']=meta
    prepared_carrier=revised_carrier.as_original().simplify(1e-6)
    actual_carrier,carrier_meta=export(prepared_carrier,DEST/'collar_guide_passage.obj',
        np.array(parent['sources']['collar:jaw_lower']['translation_source_mm']))
    sources['collar:jaw_lower']=carrier_meta
    removed=muscle-actual;added=actual-muscle
    hits=[]
    for key,row in parent['sources'].items():
        if key=='tpu:teres_major':continue
        obstacle=actual_carrier if key=='collar:jaw_lower' else load_component(row)
        v=float(A.I.intersection_volume(actual,obstacle));old=float(A.I.intersection_volume(muscle,obstacle))
        if v>.001 or old>.001:hits.append(dict(key=key,volume_mm3=v,parent_volume_mm3=old,
            new_or_increased=v>max(.001,old+.001)))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),entry_reference_sha256=sha(entry_file),
        script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        replaced_keys=['tpu:teres_major','collar:jaw_lower'],new_keys=[],all_other_206_sources_unchanged=True,
        entry_direction_source=direction.tolist(),carrier_relief_trace=trace,unresolved_outside_terminal=outside,
        dimensions=dict(terminal_radius_mm=12.,terminal_normal_mm=[-4.,.2],insertion_obstacle_inflation_mm=0.,
            flange_outer_radius_mm=10.5,central_bearing_protection_radius_mm=5.2,
            sample_step_mm=.5,fine_interval_mm=[0.,10.],fine_step_mm=.1,numerical_cleanup_mm=1e-6),
        TPU_removed_mm3=float(removed.volume()),TPU_added_mm3=float(added.volume()),
        removed_outside_distal_terminal_mm3=float((removed-zone).volume()),
        removed_within_wire_outer_2mm_guard_mm3=float((removed^protected).volume()),
        cleanup_added_mm3=float((prepared-revised).volume()),cleanup_removed_mm3=float((revised-prepared).volume()),
        carrier_removed_mm3=float((carrier-actual_carrier).volume()),
        carrier_added_mm3=float((actual_carrier-carrier).volume()),
        carrier_removed_outside_guide_zone_mm3=float(((carrier-actual_carrier)-metal_zone).volume()),
        carrier_thread_stock_lost_mm3=float(((carrier-actual_carrier)^thread_stock).volume()),
        changed_part_pairs=207,material_intersections=hits,
        no_new_or_increased_overlap=not any(row['new_or_increased'] for row in hits),
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners=parent['intended_body_owners'],adopted_as_current=False,stage_pass=False,
        whole_module_installation_verified=False,actual_drive_range_deg=None)
    save(DEST/'assembly_check.json',report)
    print({k:report[k] for k in ['TPU_removed_mm3','removed_outside_distal_terminal_mm3',
        'removed_within_wire_outer_2mm_guard_mm3','no_new_or_increased_overlap']},flush=True)
    print('Unresolved terminal extent',outside,flush=True)
    if outside or not report['no_new_or_increased_overlap'] or report['removed_within_wire_outer_2mm_guard_mm3']>.001:
        raise ValueError('Declared terminal lead-in not sufficient; see saved candidate evidence')


if __name__=='__main__':main()
