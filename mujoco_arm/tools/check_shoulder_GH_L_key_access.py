"""Finite nominal L-key envelopes in the verified bare-humerus assembly pose.

No actual screw socket, tool insertion, nut wrench or tightening validation.
Uses a conservative circular shaft envelope for a candidate 2.5 mm hex key.
"""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from open_shoulder_GH_assembly_entries import DEST
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    ordered_file=DEST/'ordered_posed_placement_checks.json';ordered=json.loads(ordered_file.read_text(encoding='utf8'))
    if ordered['assembly_sha256']!=A.hashed(source) or not ordered['all_selected_steps_clear']:
        raise ValueError('Need verified ordered placement for this source')
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    joint=json.loads((A.G1/'assembly.json').read_text(encoding='utf8'));n=np.array(joint['neck_axis']);G=np.array(joint['centre_source_mm'])
    pose=ordered['assembly_pose'];S['bone:humerus']=S['bone:humerus'].transform(A.I.rotation_about(G,pose['axis_source'],pose['angle_deg']))
    fixed_keys=set(ordered['steps'][-1]['path']['fixed_keys'])|{'GH_bolt_1'}
    x=A.D.unit(np.cross(n,[0.,-1.,0.]));y=np.cross(n,x);cases=[]
    for i,point in enumerate(joint['bolt_sites_display_mm']):
        c=np.array(point)+A.frames.SHIFT;selected=None;attempts=[]
        for short in [8.,6.,4.]:
            for angle in range(0,360,10):
                v=x*np.cos(np.radians(angle))+y*np.sin(np.radians(angle))
                start=c+n*6.75;heel=start+n*short
                body=rod(start,heel,1.5)+rod(heel,heel+v*50.,1.5)+A.I._manifold().Manifold.sphere(1.6,24).translate(heel)
                hits={k:A.I.intersection_volume(body,S[k]) for k in fixed_keys if k!=f'GH_bolt_{i}'}
                bad=[dict(part=k,volume_mm3=w) for k,w in hits.items() if w>.001]
                attempts.append(dict(short_leg_mm=short,handle_angle_deg=angle,hits_over_threshold=bad))
                if bad:continue
                turns=[]
                for direction in [-1,1]:
                    samples=[]
                    for value in np.arange(0.,21.,1.):
                        posed=body.transform(A.I.rotation_about(c,n,float(value*direction)))
                        badturn=[dict(part=k,volume_mm3=A.I.intersection_volume(posed,S[k])) for k in fixed_keys if k!=f'GH_bolt_{i}']
                        badturn=[r for r in badturn if r['volume_mm3']>.001]
                        samples.append(dict(relative_angle_deg=float(value*direction),hits_over_threshold=badturn))
                        if badturn:break
                    turns.append(dict(direction=direction,samples=samples,full_20deg_clear=len(samples)==21 and not any(r['hits_over_threshold'] for r in samples)))
                if not all(r['full_20deg_clear'] for r in turns):continue
                actual,meta=local_export(body,DEST/f'candidate_L_key_{i}.stl',c)
                if not meta['mesh_checks_pass']:raise ValueError('L key readback failed')
                selected=dict(short_leg_mm=short,long_leg_mm=50.,shaft_envelope_radius_mm=1.5,bend_envelope_radius_mm=1.6,
                              nominal_hex_across_flats_mm=2.5,handle_angle_deg=angle,source_direction=v.tolist(),
                              fixed_keys=sorted(fixed_keys-{f'GH_bolt_{i}'}),turn_checks=turns,mesh=meta)
                break
            if selected:break
        cases.append(dict(bolt=f'GH_bolt_{i}',attempts=attempts,selected=selected))
        print('L key',i,'selected',None if selected is None else (selected['short_leg_mm'],selected['handle_angle_deg']),flush=True)
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),ordered_report_sha256=A.hashed(ordered_file),
                assembly_pose=pose,cases=cases,stage_pass=False,actual_drive_range_deg=None,
                not_verified=['real screw drive recess and key engagement depth','tool insertion and operator hand clearance',
                              'nut antirotation/reaction wrench and tightening/preload','tool material, bend strength and actual torque'])
    (DEST/'L_key_envelope_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
