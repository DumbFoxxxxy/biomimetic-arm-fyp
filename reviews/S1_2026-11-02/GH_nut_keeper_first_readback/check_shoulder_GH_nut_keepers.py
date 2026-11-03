"""Independent finite pocket reaction, retention and staged keeper placement.

All nominal meshes are read back and hashed. This is not preload, strength,
thread, print tolerance or physical assembly/actuation acceptance.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_GH_nut_keepers import DEST
from check_shoulder_front_liner_keepers import dense_slide
from audit_shoulder_GH_assembly import scan


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    joint=json.loads((A.G1/'assembly.json').read_text(encoding='utf8'));G=np.array(joint['centre_source_mm']);n=np.array(joint['neck_axis'])
    definitions=data['GH_nut_keepers']['definitions'];reaction=[];capture=[]
    for definition in definitions:
        i=definition['index'];c=np.array(definition['source_origin_mm']);direction=np.array(definition['entry_direction_source'])
        for sign in [-1,1]:
            samples=[]
            for angle in np.arange(0.,10.01,.1):
                posed=S[f'GH_nut_{i}'].transform(A.I.rotation_about(c,n,float(sign*angle)))
                volume=A.I.intersection_volume(posed,S['bone:scapula'])
                samples.append(dict(angle_deg=round(float(sign*angle),6),scapula_overlap_mm3=volume))
                if volume>.001:break
            reaction.append(dict(nut=f'GH_nut_{i}',sign=sign,samples=samples,
                                 first_interference_deg=next((r['angle_deg'] for r in samples if r['scapula_overlap_mm3']>.001),None)))
        for key in [f'GH_nut_{i}',f'GH_bottom_washer_{i}']:
            fixed={'bone:scapula':S['bone:scapula'],definition['keeper_hardware']['screw']:S[definition['keeper_hardware']['screw']]}
            path=scan(S[key],direction,fixed,25.,.1)
            capture.append(dict(part=key,path=path))
    # Assembly pose is a prescribed bare-humerus condition, not an actuator.
    original_humerus=S['bone:humerus'];flex=A.frames.gh_probe_axes()['flexion']
    S['bone:humerus']=original_humerus.transform(A.I.rotation_about(G,flex,10.))
    installed={k for k in S if k.startswith(('fixture:','bench:','clamp:','liner_keeper:'))
               or k in ['bone:scapula','bone:humerus','clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    jobs=[]
    for d in definitions:
        i=d['index'];hardware=d['keeper_hardware'];v=np.array(d['keeper_axis_source'])
        jobs.extend([(f'nut_washer_{i}',[f'GH_nut_{i}',f'GH_bottom_washer_{i}'],np.array(d['entry_direction_source']),25.),
                     (f'keeper_screw_near_washer_{i}',[hardware['screw'],hardware['near_washer']],-v,15.),
                     (f'keeper_far_washer_{i}',[hardware['far_washer']],v,15.),
                     (f'keeper_nut_{i}',[hardware['nut']],v,15.)])
    jobs.extend([('cap_top_washer_0',['GH_cap_0','GH_top_washer_0'],A.D.unit(np.array([-1.,-1.,-1.])),25.),
                 ('cap_top_washer_1',['GH_cap_1','GH_top_washer_1'],n,25.),
                 ('bolt_0',['GH_bolt_0'],n,15.),('bolt_1',['GH_bolt_1'],n,15.)])
    steps=[]
    for name,keys,direction,length in jobs:
        result=dense_slide({k:S[k] for k in keys},{k:S[k] for k in sorted(installed)},direction,length)
        firsts=[dict(phase=p['name'],first_hit=next((r for r in p['samples'] if r['hits_over_threshold']),None)) for p in result['phases']]
        steps.append(dict(name=name,path=result,first_interferences=firsts,
                          not_yet_installed_keys=sorted(set(S)-installed-set(keys))))
        installed.update(keys)
        print(name,result['sampled_path_clear'],firsts if not result['sampled_path_clear'] else '',flush=True)
    fixed={k:S[k] for k in installed-{'bone:humerus'}}
    motion=[]
    for label,values in [('to_10deg',np.arange(0.,11.)),('return_neutral',np.arange(10.,-1.,-1.))]:
        samples=[]
        for angle in values:
            h=original_humerus.transform(A.I.rotation_about(G,flex,float(angle)))
            bad=[dict(part=k,volume_mm3=A.I.intersection_volume(h,s)) for k,s in fixed.items()]
            bad=[r for r in bad if r['volume_mm3']>.001];samples.append(dict(angle_deg=float(angle),hits_over_threshold=bad))
            if bad:break
        motion.append(dict(name=label,samples=samples,sampled_path_clear=len(samples)==len(values) and not any(r['hits_over_threshold'] for r in samples)))
    dirs=[];x=A.D.unit(np.cross(n,[0.,-1.,0.]));F=np.column_stack([x,np.cross(n,x),n])
    for basis in [np.eye(3),F]:
        for vec in itertools.product([-1.,0.,1.],repeat=3):
            if not any(vec):continue
            d=A.D.unit(basis@vec)
            if not any(np.linalg.norm(d-w)<1e-6 for w in dirs):dirs.append(d)
    sphere=A.I._manifold().Manifold.sphere(joint['ball_radius_mm'],40).translate(G)
    sphere_checks=[]
    for direction in dirs:
        path=scan(sphere,direction,{k:S[k] for k in ['bone:scapula','GH_cap_0','GH_cap_1']},25.,.1)
        sphere_checks.append(dict(direction_source=direction.tolist(),path=path))
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),reaction_probes=reaction,side_retention_probes=capture,
                ordered_steps=steps,bare_humerus_roundtrip=motion,finite_ball_capture=sphere_checks,
                all_source_hashes_match=True,translation_step_mm=.1,reaction_angle_step_deg=.1,humerus_angle_step_deg=1.,
                all_reaction_directions_encounter_material=all(r['first_interference_deg'] is not None for r in reaction),
                all_side_retention_directions_encounter_material=all(r['path']['first_interference'] is not None for r in capture),
                all_ordered_paths_clear=all(r['path']['sampled_path_clear'] for r in steps),
                bare_humerus_roundtrip_clear=all(r['sampled_path_clear'] for r in motion),
                all_52_ball_directions_encounter_material=all(r['path']['first_interference'] is not None for r in sphere_checks),
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['actual tightening and M2 keeper tools, thread/head geometry and preload',
                              'force/torque, nut rotation strength, keeper bending and printed fits',
                              'all escape directions, full shoulder assembly, muscle deformation and drive'])
    (DEST/'readback_keeper_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Nut rotation limits',[(r['nut'],r['first_interference_deg']) for r in reaction],flush=True)
    print('Side limits',[(r['part'],None if r['path']['first_interference'] is None else r['path']['first_interference']['distance_mm']) for r in capture],flush=True)
    print('Pose clear',report['bare_humerus_roundtrip_clear'],'ball52',report['all_52_ball_directions_encounter_material'],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
