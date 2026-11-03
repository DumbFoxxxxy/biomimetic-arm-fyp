"""Independent dense-curve, endpoint, geometry and export checks of route witnesses."""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import shoulder_free_terminal_constraints as C
import shoulder_free_strand_routes as R
import design_shoulder_tendon_link as L
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_y2_transition_stock/free_strand_routes'

def dense_metrics(record):
    control=np.array(record['control_points_source_mm']);a,b,c,d=control
    t=np.linspace(0.,1.,32769)[:,None]
    curve=(1-t)**3*a+3*t*(1-t)**2*b+3*t*t*(1-t)*c+t**3*d
    length=float(np.linalg.norm(np.diff(curve,axis=0),axis=1).sum())
    velocity=3*((1-t)**2*(b-a)+2*t*(1-t)*(c-b)+t*t*(d-c))
    accel=6*((1-t)*(c-2*b+a)+t*(d-2*c+b))
    speed=np.linalg.norm(velocity,axis=1);cross=np.linalg.norm(np.cross(velocity,accel),axis=1)
    valid=cross>1e-14;radius=float(np.min(speed[valid]**3/cross[valid])) if np.any(valid) else None
    analytic=record['minimum_cubic_radius']['radius_mm']
    if abs(length-record['cubic_arc_length_mm'])>2e-6:raise ValueError('Dense arc length mismatch')
    if analytic is not None and radius is not None and abs(radius-analytic)>max(1e-5,abs(analytic)*2e-5):raise ValueError('Dense curvature mismatch')
    return dict(dense_subdivisions=32768,dense_length_mm=length,
        quadrature_difference_mm=length-record['cubic_arc_length_mm'],dense_sampled_minimum_radius_mm=radius,
        polynomial_critical_minimum_radius_mm=analytic)

def run():
    source=DEST/'route_audit.json';audit=json.loads(source.read_text(encoding='utf8'))
    assembly_file,assembly,models,provenance=C.source_model()
    if A.hashed(assembly_file)!=audit['assembly_sha256']:raise ValueError('Changed assembly')
    if A.hashed(A.ROOT/'mujoco_arm/tools/audit_shoulder_free_strand_routes.py')!=audit['script_sha256']:raise ValueError('Changed auditor')
    for name,expected in audit['dependencies_sha256'].items():
        if A.hashed(A.ROOT/'mujoco_arm/tools'/name)!=expected:raise ValueError('Changed route dependency')
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in assembly['sources'].items()}
    curve_checks=[];endpoint_checks=[];first_reasons=[];export_checks=[];first_rigid_replays=[]
    for head in audit['heads']:
        model=next(m for m in models if m['head']==head['head'])
        cases=[('neutral',head['neutral'])]
        for direction in head['directions']:
            if not direction['reverse_candidates_match']:raise ValueError('Reverse branch failed')
            accepted=[r for r in direction['coarse_samples']+direction['boundary_samples'] if r['candidate_found']]
            cases.append((direction['pattern']+':'+str(direction['sign']),accepted[-1]))
            failed=direction['first_no_candidate']
            if failed:
                record=failed['route'];reason=record['reason']
                if reason=='exit chord exceeds assumed nominal material length':
                    if record['exit_chord_mm']<=model['nominal_active_length']+1e-6:raise ValueError('Invalid necessary-length rejection')
                elif reason=='retained straight-lead family has no positive-handle length solution':
                    if record['family_straight_lead_lower_bound_mm']<model['nominal_active_length']-1e-8:raise ValueError('Invalid route-family bound')
                elif reason=='candidate curve radius not above 0.3mm wire radius':
                    dense_metrics(record)
                if failed['failure_category']=='rigid_route_or_eye_package_interference':
                    phiA=np.deg2rad(failed['phiA_deg']);phiB=np.deg2rad(failed['phiB_deg'])
                    RA=C.rotation(model['axisA'],phiA);RB=C.rotation(model['axisB'],phiB)
                    TA=np.c_[RA,model['pinA']-RA@model['pinA']];TB=np.c_[RB,model['pinB']-RB@model['pinB']]
                    prefix=np.c_[model['prefix'],np.ones(len(model['prefix']))]@TA.T
                    suffix=np.c_[model['suffix'],np.ones(len(model['suffix']))]@TB.T
                    full=C.deduplicate(np.vstack([prefix,np.array(record['sampled_path_source_mm'])[1:-1],suffix]))
                    posed=dict(S);changed=[]
                    for end,transform in [('muscle',TA),('bone',TB)]:
                        for part in ['groove_ring','uncrimped_sleeve']:
                            key='link:'+model['head']+'/'+end+'_'+part;posed[key]=S[key].transform(transform);changed.append(key)
                    wire_key='link:'+model['head']+'/tendon_wire';posed[wire_key]=L.wire_tube(full);changed.append(wire_key)
                    hits=[];seen=set()
                    for key in changed:
                        for other,solid in posed.items():
                            if key==other or is_compliant(key) or is_compliant(other):continue
                            pair=tuple(sorted([key,other]))
                            if pair in seen:continue
                            seen.add(pair);volume=A.I.intersection_volume(posed[key],solid)
                            if volume>.001:hits.append(dict(a=key,b=other,volume_mm3=volume))
                    if not hits:raise ValueError('Independent full-route first rigid boundary disagrees')
                    first_rigid_replays.append(dict(head=model['head'],pattern=direction['pattern'],sign=direction['sign'],
                        phiA_deg=failed['phiA_deg'],phiB_deg=failed['phiB_deg'],hits_over_threshold=hits))
                first_reasons.append(dict(head=head['head'],pattern=direction['pattern'],sign=direction['sign'],
                    category=failed['failure_category'],route_reason=reason,
                    failure_is_not_a_global_motion_or_routing_infeasibility_claim=True))
        for label,state in cases:
            record=state['route'];metrics=dense_metrics(record);control=np.array(record['control_points_source_mm'])
            path=np.array(record['sampled_path_source_mm']);exitA=np.array(record['exit_A_source_mm']);exitB=np.array(record['exit_B_source_mm'])
            uA=np.array(record['exit_tangent_A']);uB=np.array(record['exit_tangent_B_toward_free_segment'])
            ta=control[1]-control[0];tb=control[2]-control[3];ta/=np.linalg.norm(ta);tb/=np.linalg.norm(tb)
            errors=[np.linalg.norm(path[0]-exitA),np.linalg.norm(path[-1]-exitB),np.linalg.norm(ta-uA),np.linalg.norm(tb-uB)]
            if max(errors)>1e-7:raise ValueError('Clamped endpoint/tangent mismatch')
            if abs(C.polyline_length(path)-model['nominal_active_length'])>1e-4:raise ValueError('Sampled route length changed')
            curve_checks.append(dict(head=head['head'],case=label,metrics=metrics))
            endpoint_checks.append(dict(head=head['head'],case=label,endpoint_and_tangent_errors=errors))
        for entry in head['diagnostic_exports']:
            scene_file=A.ROOT/entry['scene'];scene=json.loads(scene_file.read_text(encoding='utf8'))
            if scene['original_assembly_sha256']!=audit['assembly_sha256'] or len(scene['sources'])!=191:raise ValueError('Wrong diagnostic context')
            changed=[k for k,row in scene['sources'].items() if row['sha256']!=assembly['sources'][k]['sha256']]
            expected=set(['link:'+head['head']+'/'+end+'_'+part for end in ['muscle','bone'] for part in ['groove_ring','uncrimped_sleeve']]+['link:'+head['head']+'/tendon_wire'])
            if set(changed)!=expected:raise ValueError('Wrong diagnostic moving/replaced sources')
            posed=dict(S)
            for key in changed:
                row=scene['sources'][key];posed[key],_,meta=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
            rigid=[];soft=[];seen=set()
            for key in changed:
                for other,solid in posed.items():
                    if other==key:continue
                    pair=tuple(sorted([key,other]))
                    if pair in seen:continue
                    seen.add(pair);volume=A.I.intersection_volume(posed[key],solid)
                    if volume>.001:(soft if is_compliant(key) or is_compliant(other) else rigid).append(dict(a=key,b=other,volume_mm3=volume))
            if rigid or soft:raise ValueError('Exported demonstration geometry changed to an intersecting pose')
            length=C.polyline_length(np.array(entry['full_centreline_source_mm']))
            if abs(length-model['nominal_full_length'])>1e-4:raise ValueError('Full diagnostic folded length changed')
            export_checks.append(dict(head=head['head'],label=entry['label'],scene_sha256=A.hashed(scene_file),
                verified_changed_sources=changed,parts_in_context=191,changed_vs_context_pairs=len(seen),
                readback_rigid_hits=rigid,readback_undeformed_TPU_hits=soft,full_wire_polyline_length_error_mm=length-model['nominal_full_length']))
    report=dict(scope=__doc__,audit_sha256=A.hashed(source),assembly_sha256=A.hashed(assembly_file),
        script_sha256=A.hashed(Path(__file__).resolve()),dense_curve_checks=curve_checks,
        clamped_endpoint_checks=endpoint_checks,first_failure_classifications=first_reasons,diagnostic_export_checks=export_checks,
        independent_first_rigid_boundary_replays=first_rigid_replays,
        all_dense_curves_and_exports_match=True,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['global shortest route or all curve families','installed length and bending stiffness/allowable radius','TPU/cap/humerus equilibrium and driven motion'])
    (DEST/'verification_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Dense curve cases',len(curve_checks),'; independent readback scenes',len(export_checks),'all match',flush=True)

if __name__=='__main__':run()
