"""Adopt the checked geometry research candidate, retaining all prior inputs."""
import copy
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'
OLD=D/'s1_y_muscle_pair';NEW=D/'s1_y2_transition_stock'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

pointer_file=ROOT/'design/S1_CURRENT_CANDIDATE.json'
before=Path(__file__).resolve().parent/'before/design/S1_CURRENT_CANDIDATE.json'
if sha(pointer_file)!=sha(before):raise ValueError('Selection changed since review snapshot')
selection=read(pointer_file);assembly=read(NEW/'assembly_check.json')
if selection['assembly_sha256']!=assembly['parent_sha256']:raise ValueError('Wrong selected parent')
if set(assembly['changed_keys'])!={'tpu:infraspinatus','tpu:subscapularis'}:raise ValueError('Unexpected changes')
if not all(sha(ROOT/r['file'])==r['sha256'] for r in assembly['sources'].values()):raise ValueError('Changed source')
if not assembly['no_rigid_overlap_over_threshold'] or not assembly['no_wire_overlap_over_threshold']:raise ValueError('Static fit failed')
for rel in ['bare_motion_review/bare_joint_motion_checks.json','readback_keeper_checks.json']:
    if read(NEW/rel)['assembly_sha256']!=sha(NEW/'assembly_check.json'):raise ValueError('Stale new report')
keeper=read(NEW/'readback_keeper_checks.json')
if not all(keeper[k] for k in ['all_source_hashes_match','all_ordered_paths_clear','bare_humerus_roundtrip_clear','all_52_ball_directions_encounter_material']):raise ValueError('Keeper checks failed')
if not all(m['verification']['complete_pocket_resolution_geometry_pass'] for m in assembly['muscles']):raise ValueError('Muscle checks failed')

wire=copy.deepcopy(read(OLD/'wire_and_anchor_plan.json'))
wire['scope']='S1-Y2 nominal wire/interface references, unchanged paths; no validated drive stroke or load.'
wire['assembly_sha256']=sha(NEW/'assembly_check.json')
for m in wire['muscles']:m['TPU_mesh']=assembly['sources']['tpu:'+m['module']]
write(NEW/'wire_and_anchor_plan.json',wire)
template=copy.deepcopy(read(OLD/'single_muscle_observation_template.json'))
if any(r['measured_test_completed'] or r['observations'] for r in template['specimens']):raise ValueError('Do not relabel measured records')
for r in template['specimens']:r['source_mesh_sha256']=assembly['sources']['tpu:'+r['module']]['sha256']
write(NEW/'single_muscle_observation_template.json',template)

selection['candidate']='S1-Y2: compact GH skeleton with two conforming Y bellies and explicitly infilled closed transition pockets'
selection['parent_research_report']=selection['assembly_report'];selection['parent_research_sha256']=selection['assembly_sha256']
for file_key,hash_key,rel in [('assembly_report','assembly_sha256','assembly_check.json'),('bare_motion_report','bare_motion_report_sha256','bare_motion_review/bare_joint_motion_checks.json'),('keeper_report','keeper_report_sha256','readback_keeper_checks.json')]:
    p=NEW/rel;selection[file_key]=p.relative_to(ROOT).as_posix();selection[hash_key]=sha(p)
for name,rel in [('transition_stock_generation_report','transition_trials.json'),('transition_stock_independent_checks','independent_checks.json')]:
    p=OLD/'transition_stock_comparison'/rel;selection[name]=p.relative_to(ROOT).as_posix();selection[name+'_sha256']=sha(p)
selection['limits'].append('explicit closed-pocket infill changes local stiffness/folding; mechanical and manufacturing effects remain unmeasured')
selection['stage_pass']=False;selection['actual_drive_range_deg']=None
write(pointer_file,selection)
print('S1-Y2 selected as checked geometry research candidate; not manufacturing or drive release.')
