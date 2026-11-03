"""Ordered bare-shoulder GH subassembly placement after local entry revisions.

Checks each step against all previously installed nominal objects and performs
independent forward/inverse samples. Manual temporary holding, nut reaction,
thread tightening and full-muscle installation remain separate requirements.
"""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from open_shoulder_GH_assembly_entries import DEST
from check_shoulder_front_liner_keepers import dense_slide
from build_shoulder_shaft_collar import rod


def run():
    source=DEST/'assembly_check.json';data=json.loads(source.read_text(encoding='utf8'))
    solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in data['sources'].items()}
    jointfile=A.G1/'assembly.json';joint=json.loads(jointfile.read_text(encoding='utf8'))
    n=np.array(joint['neck_axis'])
    installed={k for k in solids if k.startswith(('fixture:','bench:','clamp:','liner_keeper:'))
               or k in ['bone:scapula','bone:humerus','clavicle:fixed_support'] or k.startswith('clavicle:M4_')}
    steps=[]
    entries=data['GH_assembly_entries']['definitions']
    jobs=[(r['name'],r['moving_keys'],np.array(r['direction_source']),r['maximum_translation_mm']) for r in entries[:2]]
    jobs += [('cap_and_top_washer_0',['GH_cap_0','GH_top_washer_0'],A.D.unit(np.array([-1.,-1.,-1.])),25.),
             ('cap_and_top_washer_1',['GH_cap_1','GH_top_washer_1'],n,25.),
             ('bolt_0',['GH_bolt_0'],n,15.),('bolt_1',['GH_bolt_1'],n,15.)]
    for name,keys,direction,maximum in jobs:
        result=dense_slide({k:solids[k] for k in keys},{k:solids[k] for k in sorted(installed)},direction,maximum)
        firsts=[dict(phase=p['name'],first_hit=next((s for s in p['samples'] if s['hits_over_threshold']),None)) for p in result['phases']]
        steps.append(dict(name=name,placement_order=len(steps)+1,path=result,first_interferences=firsts,
                          not_yet_installed_keys=sorted(set(solids)-installed-set(keys))))
        installed.update(keys)
        print(name,result['sampled_path_clear'],firsts if not result['sampled_path_clear'] else '',flush=True)
    tools=[]
    for i,point in enumerate(joint['bolt_sites_display_mm']):
        c=np.array(point)+A.frames.SHIFT
        shank=rod(c+n*6.75,c+n*36.75,1.5)
        handle=rod(c+n*36.75,c+n*56.75,8.)
        envelope=shank+handle
        hits={k:A.I.intersection_volume(envelope,solids[k]) for k in installed if k!=f'GH_bolt_{i}'}
        tools.append(dict(bolt=f'GH_bolt_{i}',scope='one final-pose hypothetical cylindrical tool, not a real tool or swept insertion',
                          shank_radius_mm=1.5,shank_length_mm=30.,handle_radius_mm=8.,handle_length_mm=20.,
                          fixed_keys=sorted(installed-{f'GH_bolt_{i}'}),
                          hits_over_threshold=[dict(part=k,volume_mm3=v) for k,v in hits.items() if v>.001]))
    report=dict(scope=__doc__,assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),joint_reference_sha256=A.hashed(jointfile),
                dense_slide_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/check_shoulder_front_liner_keepers.py'),
                steps=steps,selected_tool_envelopes=tools,all_source_hashes_match=True,step_mm=.1,
                all_selected_steps_clear=all(s['path']['sampled_path_clear'] for s in steps),
                all_selected_tool_envelopes_clear=not any(t['hits_over_threshold'] for t in tools),
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['temporary holding and retention/anti-rotation of the side-entered nuts and washers',
                              'actual driver/head engagement, wrench access, tool insertion and preload',
                              'scapula clamping after bone revisions, full skeleton/muscle installation',
                              'multi-direction capture, printed clearance, strength and muscle drive'])
    (DEST/'ordered_placement_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Hypothetical tool hits',[t['hits_over_threshold'] for t in tools],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
