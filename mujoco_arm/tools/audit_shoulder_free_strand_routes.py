"""Full folded-wire route witnesses for prescribed eye rotations only.

Caps, TPU, humerus and fixtures stay neutral. Both real eye packages may
rotate. The entire own tendon wire is replaced with length-matched geometry;
no central strand is omitted. Nominal material length is an explicit
assumption, not an installed measurement, taut equilibrium or drive range.
"""
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
from build_shoulder_Y_print_inspection import write_readback

DEST=A.DESIGN/'s1_y2_transition_stock/free_strand_routes'
PATTERNS={'cap_eye_only':(1,0),'bone_eye_only':(0,1),'same_signed_eyes':(1,1),'opposite_signed_eyes':(1,-1)}

def nonlocal_self_distance(path,cap=1.):
    p=np.asarray(path);lo=np.minimum(p[:-1],p[1:]);hi=np.maximum(p[:-1],p[1:])
    arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))];best=cap;pair=None
    for i in range(len(p)-1):
        eligible=(arc[:-1]-arc[i+1]>1.2)
        eligible&=np.all(lo<=hi[i]+best,axis=1)&np.all(hi>=lo[i]-best,axis=1)
        for j in np.flatnonzero(eligible):
            distance=A.D.segment_distance(p[i],p[i+1],p[j],p[j+1])
            if distance<best:best=distance;pair=[int(i),int(j)]
    return dict(distance_capped_1mm=float(best),segment_pair=pair,excluded_arclength_neighbours_mm=1.2,
        self_contact_over_threshold=best<.6-1e-6)

def run():
    assembly_file,assembly,heads,provenance=C.source_model()
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in assembly['sources'].items()}
    DEST.mkdir(parents=True,exist_ok=True);head_reports=[]
    for model in heads:
        head=model['head'];wire_key='link:'+head+'/tendon_wire'
        movable=['link:'+head+'/'+end+'_'+part for end in ['muscle','bone'] for part in ['groove_ring','uncrimped_sleeve']]
        def evaluate(a,b):
            phiA=np.deg2rad(a);phiB=np.deg2rad(b)
            path,route=R.route(model,phiA,phiB)
            row=dict(phiA_deg=round(float(a),6),phiB_deg=round(float(b),6),route=route,
                candidate_found=False,tension_N=None,pose_is_not_solved_equilibrium=True)
            if path is None:row['failure_category']='necessary_length_bound' if route['conditional_impossibility_for_prescribed_pose'] else 'declared_route_family';return row,None
            TA=A.I.rotation_about(model['pinA'],model['axisA'],a);TB=A.I.rotation_about(model['pinB'],model['axisB'],b)
            prefix=np.c_[model['prefix'],np.ones(len(model['prefix']))]@TA.T
            suffix=np.c_[model['suffix'],np.ones(len(model['suffix']))]@TB.T
            full=C.deduplicate(np.vstack([prefix,path[1:-1],suffix]))
            row['full_wire_length_error_mm']=C.polyline_length(full)-model['nominal_full_length']
            row['nonlocal_self_clearance']=nonlocal_self_distance(full)
            if row['nonlocal_self_clearance']['self_contact_over_threshold']:
                row['failure_category']='nonlocal_wire_self_contact';return row,None
            newwire=L.wire_tube(full);row['wire_manifold_status']=str(newwire.status())
            if newwire.is_empty() or 'NoError' not in str(newwire.status()):
                row['failure_category']='wire_sweep_mesh';return row,None
            posed=dict(S)
            for key in movable:posed[key]=S[key].transform(TA if '/muscle_' in key else TB)
            posed[wire_key]=newwire;changed=movable+[wire_key];seen=set();rigid=[];soft=[]
            for key in changed:
                for other,solid in posed.items():
                    if key==other:continue
                    pair=tuple(sorted([key,other]))
                    if pair in seen:continue
                    seen.add(pair);volume=A.I.intersection_volume(posed[key],solid)
                    if volume>.001:(soft if is_compliant(key) or is_compliant(other) else rigid).append(dict(a=key,b=other,volume_mm3=volume))
            row.update(rigid_hits_over_threshold=rigid,undeformed_TPU_hits_over_threshold=soft,
                changed_vs_context_pairs_checked=len(seen),candidate_found=not rigid,
                failure_category=None if not rigid else 'rigid_route_or_eye_package_interference')
            return row,dict(full=full,posed={k:posed[k] for k in changed})
        neutral,neutral_geom=evaluate(0.,0.)
        if not neutral['candidate_found']:raise ValueError('Neutral clamped route family failed: '+head+' '+str(neutral))
        directions=[];witnesses=[(neutral,neutral_geom)]
        for label,(ca,cb) in PATTERNS.items():
            for sign in [1,-1]:
                clear=[];coarse=[];fine=[];first=None
                for magnitude in range(181):
                    alpha=sign*magnitude;row,geometry=evaluate(ca*alpha,cb*alpha);row['sweep_parameter_deg']=float(alpha);coarse.append(row)
                    if not row['candidate_found']:first=row;break
                    clear.append(row)
                    if not row.get('undeformed_TPU_hits_over_threshold'):witnesses.append((row,geometry))
                if first is not None and clear:
                    base=clear[-1]['sweep_parameter_deg']
                    for i in range(1,11):
                        alpha=base+sign*i*.1;row,geometry=evaluate(ca*alpha,cb*alpha);row['sweep_parameter_deg']=round(float(alpha),6);fine.append(row)
                        if not row['candidate_found']:first=row;break
                        clear.append(row)
                reverse=[]
                for r in reversed(clear):
                    check,_=evaluate(r['phiA_deg'],r['phiB_deg'])
                    reverse.append(dict(phiA_deg=r['phiA_deg'],phiB_deg=r['phiB_deg'],candidate_found=check['candidate_found'],failure_category=check['failure_category']))
                row=dict(pattern=label,coefficients=[ca,cb],sign=sign,
                    last_candidate_parameter_deg=clear[-1]['sweep_parameter_deg'] if clear else None,
                    first_no_candidate=first,coarse_samples=coarse,boundary_samples=fine,
                    independently_rebuilt_reverse=reverse,reverse_candidates_match=all(r['candidate_found'] for r in reverse),
                    not_a_global_feasible_angle_range=True)
                directions.append(row)
                print(head,label,sign,'last candidate',row['last_candidate_parameter_deg'],'first category',first['failure_category'] if first else 'scan limit',flush=True)
        # Export the neutral and an interior demonstration witness without
        # undeformed TPU intersections. These are prescribed diagnostic poses.
        witnesses.sort(key=lambda item:max(abs(item[0]['phiA_deg']),abs(item[0]['phiB_deg'])))
        selected=[('neutral',neutral,neutral_geom)]
        if len(witnesses)>1:
            row,geometry=witnesses[len(witnesses)//2];selected.append(('prescribed_eye_witness',row,geometry))
        exports=[]
        for name,state,geometry in selected:
            folder=DEST/head/name;folder.mkdir(parents=True,exist_ok=True);sources=dict(assembly['sources']);meta_list={}
            for key,solid in geometry['posed'].items():
                origin=np.array(assembly['sources'][key]['translation_source_mm']);local,meta=write_readback(solid.translate(-origin),folder/(key.split('/')[-1]+'.stl'))
                meta['translation_source_mm']=origin.tolist();sources[key]=meta;meta_list[key]=meta
                if meta['boundary_edges'] or meta['nonmanifold_edges'] or meta['positive_components']!=1:raise ValueError('Diagnostic export invalid')
            scene=dict(scope=__doc__,sources=sources,n_instances=len(sources),phiA_deg=state['phiA_deg'],phiB_deg=state['phiB_deg'],
                original_assembly_sha256=A.hashed(assembly_file),central_strand_present=True,prescribed_pose_not_equilibrium=True)
            (folder/'scene.json').write_text(json.dumps(scene,indent=2)+'\n',encoding='utf8')
            exports.append(dict(label=name,scene=(folder/'scene.json').relative_to(A.ROOT).as_posix(),phiA_deg=state['phiA_deg'],phiB_deg=state['phiB_deg'],
                full_centreline_source_mm=geometry['full'].tolist(),exported_sources=meta_list))
        chord=(model['exitB']-model['exitA']);unit=chord/np.linalg.norm(chord)
        head_reports.append(dict(head=head,nominal_free_strand_length_assumption_mm=model['nominal_active_length'],
            exit_axis_offset_mm=float(chord@model['axisA']),neutral_chord_out_of_eye_plane_deg=float(np.degrees(np.arcsin(np.clip(abs(unit@model['axisA']),0,1)))),
            neutral=neutral,directions=directions,diagnostic_exports=exports))
    report=dict(scope=__doc__,assembly_source=assembly_file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        source_provenance=provenance,script_sha256=A.hashed(Path(__file__).resolve()),
        dependencies_sha256={name:A.hashed(A.ROOT/'mujoco_arm/tools'/name) for name in ['shoulder_free_terminal_constraints.py','shoulder_free_strand_routes.py','build_shoulder_Y_print_inspection.py','design_shoulder_tendon_link.py']},
        heads=head_reports,angular_step_deg=1.,boundary_step_deg=.1,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['installed length, no-slip grip, bending/friction and preload','shortest route or taut/slack equilibrium','TPU/cap/humerus coupled motion','physical bend radius, fatigue or manufacturing'])
    (DEST/'route_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
