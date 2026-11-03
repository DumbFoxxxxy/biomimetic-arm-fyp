"""Finite-difference and independent boundary replay of free-end diagnostics."""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import shoulder_free_terminal_constraints as C
import design_shoulder_tendon_link as L
from check_shoulder_front_liner_keepers import is_compliant
from build_shoulder_shaft_collar import rod

DEST=A.DESIGN/'s1_y2_transition_stock/free_terminal_kinematics'

def run():
    audit_file=DEST/'kinematics_audit.json';audit=json.loads(audit_file.read_text(encoding='utf8'))
    assembly_file,assembly,models,provenance=C.source_model()
    if audit['assembly_sha256']!=A.hashed(assembly_file):raise ValueError('Changed assembly')
    if audit['constraint_helper_sha256']!=A.hashed(Path(C.__file__).resolve()):raise ValueError('Changed constraints')
    if audit['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/audit_shoulder_free_terminal_kinematics.py'):raise ValueError('Changed auditor')
    results=[];input_guards=[]
    example=models[0]
    cases=[('unknown installed length',None,'measured',None,True),
           ('unverified straight route',example['nominal_active_length'],'explicit_assumption',None,False),
           ('folded path not declared as measurement',example['nominal_full_length'],'nominal',None,True),
           ('route shorter than chord',example['nominal_active_length'],'explicit_assumption',1.,False)]
    for label,installed,provenance_name,route_length,clear in cases:
        try:C.evaluate_installed_strand(example,installed,length_provenance=provenance_name,
            clear_routed_length_mm=route_length,straight_chord_clear_verified=clear)
        except ValueError as e:input_guards.append(dict(case=label,expected_rejection=True,message=str(e)))
        else:raise ValueError('Unsafe length-state input accepted: '+label)
    diagnostic=C.evaluate_installed_strand(example,example['nominal_active_length'],
        length_provenance='explicit_assumption',straight_chord_clear_verified=True)
    if diagnostic['length_provenance']!='explicit_assumption' or diagnostic['tension_N'] is not None or diagnostic['stage_pass']:
        raise ValueError('Diagnostic claimed measured force or stage pass')
    for model in models:
        for index in range(8):
            x=model['pinA']+np.array([.2*index,-.13*index,.07*index])
            Rc=C.rotation([1.,2.,-1.],.04*index);Rh=C.rotation([.5,-.8,1.],-.07*index)
            pa=.06*index;pb=-.035*index
            length,J=C.chord_and_jacobian(model,x,Rc,Rh,pa,pb)
            numerical=[]
            for k in range(8):
                step=1e-4 if k<3 else 1e-6;values=[]
                for sign in [1.,-1.]:
                    xp=x.copy();rp=Rc.copy();a=pa;b=pb
                    if k<3:xp[k]+=sign*step
                    elif k<6:rp=C.rotation(np.eye(3)[k-3],sign*step)@Rc
                    elif k==6:a+=sign*step
                    else:b+=sign*step
                    qa,qb,_=C.endpoints(model,xp,rp,Rh,a,b)
                    values.append(float(np.sqrt(np.dot(qa-qb,qa-qb))))
                numerical.append((values[0]-values[1])/(2*step))
            error=float(np.max(abs(J-numerical)))
            if error>2e-6:raise ValueError('Jacobian finite difference mismatch')
            results.append(dict(head=model['head'],synthetic_algebra_case=index,max_derivative_error=error,
                not_a_collision_or_equilibrium_verified_pose=True))
    solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in assembly['sources'].items()}
    boundaries=[];soft_onsets=[]
    for row in audit['local_swivels']:
        head=row['head'];label='muscle' if row['end']=='cap' else 'bone';pin=np.array(row['pin_centre_source_mm']);axis=np.array(row['axis_source'])
        ring='link:'+head+'/'+label+'_groove_ring';sleeve='link:'+head+'/'+label+'_uncrimped_sleeve'
        moving=solids[ring]+solids[sleeve]+L.wire_tube(np.array(row['local_wire_centreline_source_mm']))
        excluded=set([ring,sleeve,row['absent_entire_source_key']])
        fixed={k:s for k,s in solids.items() if k not in excluded and not is_compliant(k)}
        for direction in row['directions']:
            angles=[('last_clear',direction['last_clear_angle_deg'])]
            if direction['first_rigid_interference'] is not None:angles.append(('first_hit',direction['first_rigid_interference']['angle_deg']))
            for kind,angle in angles:
                if angle is None:raise ValueError('Zero pose invalid')
                # Use a separately formed homogeneous pose about the recorded
                # pin rather than the scanner's rotation_about implementation.
                R=C.rotation(axis,np.deg2rad(angle));posed=moving.transform(np.c_[R,pin-R@pin])
                hits=[]
                for key,s in fixed.items():
                    volume=A.I.intersection_volume(posed,s)
                    if volume>.001:hits.append(dict(part=key,volume_mm3=volume))
                consistent=(not hits) if kind=='last_clear' else bool(hits)
                if not consistent:raise ValueError('Independent boundary replay disagrees')
                boundaries.append(dict(head=head,end=row['end'],kind=kind,angle_deg=angle,hits_over_threshold=hits,consistent=True))
            first_soft=next((r for r in direction['coarse_samples'] if r['undeformed_TPU_hits_over_threshold']),None)
            if first_soft is not None:
                sign=direction['sign'];base=first_soft['angle_deg']-sign;samples=[]
                for k in range(1,11):
                    angle=base+sign*k*.1;R=C.rotation(axis,np.deg2rad(angle));posed=moving.transform(np.c_[R,pin-R@pin]);hits=[]
                    for key,s in solids.items():
                        if key in excluded or not is_compliant(key):continue
                        volume=A.I.intersection_volume(posed,s)
                        if volume>.001:hits.append(dict(part=key,volume_mm3=volume))
                    samples.append(dict(angle_deg=round(float(angle),6),undeformed_TPU_hits_over_threshold=hits))
                    if hits:break
                soft_onsets.append(dict(head=head,end=row['end'],sign=sign,step_deg=.1,refined_samples=samples,
                    interpretation='undeformed contact onset; neither flexible failure nor flexible pass'))
    chord_probes=[]
    for model in models:
        a,b,_=C.endpoints(model);wire=rod(a,b,.3);hits=[]
        for key,s in solids.items():
            if key=='link:'+model['head']+'/tendon_wire':continue
            volume=A.I.intersection_volume(wire,s)
            if volume>.001:hits.append(dict(part=key,volume_mm3=volume,compliant=is_compliant(key)))
        chord_probes.append(dict(head=model['head'],wire_diameter_mm=.6,parts_checked=len(solids)-1,
            hits_over_threshold=hits,neutral_straight_chord_has_no_rigid_hit=not any(not h['compliant'] for h in hits),
            interpretation='neutral straight-strand probe only, not variable-pose route solving or preload evidence'))
    for m in audit['models']:
        if m['local_path_partition_length_error_mm']>1e-8:raise ValueError('Wire partition does not conserve path length')
        if m['taut_chord_constraint_rank']!=1 or m['remaining_kinematic_unknowns_if_taut']!=7:raise ValueError('Wrong constraint count')
    report=dict(scope=__doc__,audit_sha256=A.hashed(audit_file),assembly_sha256=A.hashed(assembly_file),
        script_sha256=A.hashed(Path(__file__).resolve()),derivative_cases=results,boundary_replays=boundaries,
        refined_undeformed_TPU_contact_onsets=soft_onsets,neutral_chord_wire_probes=chord_probes,
        conditional_length_input_guards=input_guards,explicit_assumption_example=diagnostic,
        all_derivative_cases_pass=True,all_independent_boundary_replays_match=True,
        no_cap_or_wire_pose_is_a_solved_equilibrium=True,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'verification_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Finite derivatives',len(results),'max error',max(r['max_derivative_error'] for r in results),
          '; independent end-boundary replays',len(boundaries),'all match',flush=True)

if __name__=='__main__':run()
