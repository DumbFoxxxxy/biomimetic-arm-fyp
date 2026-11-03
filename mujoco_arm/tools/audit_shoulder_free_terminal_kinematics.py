"""Neutral end-package swivel clearance and source-based tendon constraints.

Only one eye/crimp/local wire package rotates at a time; the cap, collar,
other endpoint and undeformed context stay at their neutral poses. The own
central tendon strand is absent because its changed route is not solved.
These are staged local checks, never integrated muscle-driven motion.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import design_shoulder_tendon_link as L
import shoulder_free_terminal_constraints as C
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_y2_transition_stock/free_terminal_kinematics'

def scan_package(head,end,S):
    label='muscle' if end=='cap' else 'bone';pin=head['pinA'] if end=='cap' else head['pinB']
    axis=head['axisA'] if end=='cap' else head['axisB'];path=head['prefix'] if end=='cap' else head['suffix']
    ring='link:'+head['head']+'/'+label+'_groove_ring';sleeve='link:'+head['head']+'/'+label+'_uncrimped_sleeve'
    tendon='link:'+head['head']+'/tendon_wire'
    components=[S[ring],S[sleeve],L.wire_tube(path)]
    moving=A.I._manifold().Manifold.batch_boolean(components,A.I._manifold().OpType.Add)
    ignored=[ring,sleeve,tendon]
    obstacles={k:v for k,v in S.items() if k not in ignored}
    bounds={k:np.array(v.bounding_box()) for k,v in obstacles.items()}
    def evaluate(angle):
        posed=moving.transform(A.I.rotation_about(pin,axis,angle));bb=np.array(posed.bounding_box());hits=[];soft=[]
        for key,body in obstacles.items():
            ob=bounds[key]
            if np.any(bb[:3]>ob[3:]) or np.any(ob[:3]>bb[3:]):continue
            volume=A.I.intersection_volume(posed,body)
            if volume>.001:(soft if is_compliant(key) else hits).append(dict(part=key,volume_mm3=volume))
        return dict(angle_deg=round(float(angle),6),rigid_hits_over_threshold=hits,undeformed_TPU_hits_over_threshold=soft)
    zero=evaluate(0.);directions=[]
    for sign in [1,-1]:
        coarse=[];fine=[];clear=[];first=None
        if not zero['rigid_hits_over_threshold']:
            for magnitude in range(181):
                row=evaluate(sign*magnitude);coarse.append(row)
                if row['rigid_hits_over_threshold']:first=row;break
                clear.append(row)
            if first is not None and clear:
                previous=clear[-1]['angle_deg']
                for i in range(1,11):
                    row=evaluate(previous+sign*i*.1);fine.append(row)
                    if row['rigid_hits_over_threshold']:first=row;break
                    clear.append(row)
        else:first=zero
        reverse=[evaluate(r['angle_deg']) for r in reversed(clear)]
        last=clear[-1]['angle_deg'] if clear else None
        gap=None
        if last is not None:
            posed=moving.transform(A.I.rotation_about(pin,axis,last));pb=np.array(posed.bounding_box());minimum=1.;nearest=[]
            for key,body in obstacles.items():
                if is_compliant(key):continue
                bb=bounds[key];delta=np.maximum(np.maximum(pb[:3]-bb[3:],bb[:3]-pb[3:]),0.)
                if np.linalg.norm(delta)>minimum:continue
                value=A.I.gap(posed,body,minimum)
                if value<=minimum:minimum=value;nearest=[key]
            gap=dict(gap_capped_1mm=float(minimum),nearest_parts=nearest)
        directions.append(dict(sign=sign,last_clear_angle_deg=last,first_rigid_interference=first,
            reached_scan_limit=first is None,coarse_samples=coarse,boundary_samples=fine,independent_reverse_samples=reverse,
            reverse_rigid_clear=all(not r['rigid_hits_over_threshold'] for r in reverse),last_clear_minimum_gap=gap))
        print(head['head'],label,sign,'last clear',last,'first',first['rigid_hits_over_threshold'] if first else None,flush=True)
    return dict(head=head['head'],end=end,pin_centre_source_mm=pin.tolist(),axis_source=axis.tolist(),
        moving_actual_source_keys=[ring,sleeve],local_wire_centreline_source_mm=path.tolist(),
        absent_entire_source_key=tendon,absent_role='own central strand and own opposite local wire: pose-dependent routing unsolved',
        rigid_obstacles_checked=sum(not is_compliant(k) for k in obstacles),soft_obstacles_checked=sum(is_compliant(k) for k in obstacles),
        zero_pose=zero,directions=directions,angular_step_deg=1.,boundary_step_deg=.1,
        physical_pin_revolute_and_no_slip_crimp_assumptions_unverified=True,stage_pass=False)

def run():
    assembly_file,assembly,heads,provenance=C.source_model()
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in assembly['sources'].items()}
    models=[];swivels=[]
    for h in heads:
        chord,J=C.chord_and_jacobian(h)
        models.append(dict(head=h['head'],free_cap_key='free_terminal:'+h['head'],
            neutral_cap_pin_source_mm=h['pinA'].tolist(),neutral_bone_pin_source_mm=h['pinB'].tolist(),
            cap_pin_axis_source=h['axisA'].tolist(),bone_pin_axis_source=h['axisB'].tolist(),
            neutral_crimp_outer_exit_A_source_mm=h['exitA'].tolist(),neutral_crimp_outer_exit_B_source_mm=h['exitB'].tolist(),
            full_folded_wire_length_mm=h['nominal_full_length'],nominal_free_strand_material_length_mm=h['nominal_active_length'],
            neutral_exit_chord_mm=chord,nominal_route_length_minus_chord_mm=h['nominal_active_length']-chord,
            route_surplus_is_NOT_measured_slack_or_preload=True,active_strand_centreline_source_mm=h['active_path'].tolist(),
            local_path_partition_length_error_mm=abs(C.polyline_length(h['prefix'])+C.polyline_length(h['suffix'])+h['nominal_active_length']-h['nominal_full_length']),
            exit_location_errors_mm=h['exit_location_errors_mm'],exit_polyline_parameters=h['exit_polyline_parameters'],
            verified_hardware_keys=h['verified_hardware_keys'],
            unknowns_for_prescribed_humerus=['free-cap translation 3','free-cap orientation 3','cap eye swivel 1','bone eye swivel 1'],
            chord_length_jacobian_order=['cap_x_mm','cap_y_mm','cap_z_mm','cap_rx_rad','cap_ry_rad','cap_rz_rad','cap_swivel_rad','bone_swivel_rad'],
            neutral_chord_jacobian=J.tolist(),taut_chord_constraint_rank=int(np.linalg.matrix_rank(J.reshape(1,-1))),
            remaining_kinematic_unknowns_if_taut=8-int(np.linalg.matrix_rank(J.reshape(1,-1))),
            slack_chord_constraint_rank=0,condition='g=L_installed-length(clear routed strand)>=0; tension>=0; tension*g=0',
            installed_free_strand_length_mm=None,crimp_no_slip_verified=False,physical_preload_N=None,
            unresolved=['deformed TPU cap contact and restoring wrench','central strand route/obstacle contact and eye/crimp friction','installed lengths, crimp grip and external loads']))
        for end in ['cap','bone']:swivels.append(scan_package(h,end,S))
    DEST.mkdir(parents=True,exist_ok=True)
    report=dict(scope=__doc__,assembly_source=assembly_file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        source_provenance=provenance,script_sha256=A.hashed(Path(__file__).resolve()),
        constraint_helper_sha256=A.hashed(Path(C.__file__).resolve()),models=models,local_swivels=swivels,
        actual_drive_range_deg=None,stage_pass=False,
        not_verified=['complete cap equilibrium or driven arm range','central tendon routing during any swivel','TPU deformation, crimp grip, physical preload, friction, printing or strength'])
    (DEST/'kinematics_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
