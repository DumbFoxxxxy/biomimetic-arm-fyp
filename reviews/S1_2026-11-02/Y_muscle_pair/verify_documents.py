"""Read-only source/document consistency checks; writes its verification summary."""
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / '02_人体骨骼晶格肌肉方案/输出/shoulder_design'
Y = D / 's1_y_muscle_pair'
HERE = Path(__file__).resolve().parent
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    return json.loads(p.read_text(encoding='utf-8'))

selection = read(ROOT / 'design/S1_CURRENT_CANDIDATE.json')
active_report = ROOT / selection['assembly_report']
ACTIVE = active_report.parent
report = read(active_report)
parent = read(ROOT / selection['bone_reference_report'])
generation = read(Y / 'candidate_generation.json')
errors = []
for file_key, hash_key in [('assembly_report', 'assembly_sha256'), ('bare_motion_report', 'bare_motion_report_sha256'), ('keeper_report', 'keeper_report_sha256'), ('bone_reference_report', 'bone_reference_sha256'), ('Y_pair_generation_report', 'Y_pair_generation_sha256')]:
    if sha(ROOT / selection[file_key]) != selection[hash_key]:
        errors.append('Selected report hash: ' + file_key)
for key, source in report['sources'].items():
    if sha(ROOT / source['file']) != source['sha256']:
        errors.append('Source hash: ' + key)
changed = [key for key in report['sources'] if report['sources'][key]['sha256'] != parent['sources'][key]['sha256'] or report['sources'][key].get('translation_source_mm') != parent['sources'][key].get('translation_source_mm')]
if set(changed) != {'tpu:infraspinatus', 'tpu:subscapularis'}:
    errors.append('Unexpected changed source keys')
for rel, expected in generation['generation_dependency_sha256'].items():
    if sha(ROOT / rel) != expected:
        errors.append('Generation dependency: ' + rel)
active_script = 'mujoco_arm/tools/build_shoulder_Y2_transition_assembly.py' if ACTIVE.name == 's1_y2_transition_stock' else 'mujoco_arm/tools/check_shoulder_y_muscle_pair.py'
for rel, expected in [('mujoco_arm/tools/build_shoulder_y_muscle_pair.py', generation['script_sha256']), (active_script, report['script_sha256'])]:
    if sha(ROOT / rel) != expected:
        errors.append('Script hash: ' + rel)

for filename in ['bare_motion_review/bare_joint_motion_checks.json', 'readback_keeper_checks.json', 'wire_and_anchor_plan.json']:
    if read(ACTIVE / filename)['assembly_sha256'] != sha(active_report):
        errors.append('Stale assembly reference: ' + filename)
wire = read(ACTIVE / 'wire_and_anchor_plan.json')
if sha(ROOT / wire['master_routes']) != wire['master_routes_sha256']:
    errors.append('Stale master route')
for specimen in read(ACTIVE / 'single_muscle_observation_template.json')['specimens']:
    if specimen['source_mesh_sha256'] != report['sources']['tpu:' + specimen['module']]['sha256']:
        errors.append('Stale specimen: ' + specimen['module'])
    if specimen['measured_test_completed']:
        errors.append('Blank template claims measurement')

before_brief = (HERE / 'before/design/DESIGN_BRIEF.md').read_bytes()
now_brief = (ROOT / 'design/DESIGN_BRIEF.md').read_bytes()
marker = '## 历史记录阅读说明'.encode('utf-8')
history_preserved = before_brief.split(marker, 1)[1] == now_brief.split(marker, 1)[1]
if not history_preserved:
    errors.append('Historical task brief changed')
baseline = ROOT / 'design/CURRENT_BASELINE.md'
old_baseline = (HERE / 'before/design/CURRENT_BASELINE.md').read_text(encoding='utf-8')
new_baseline = baseline.read_text(encoding='utf-8')
indicator_rows = lambda s: [line for line in s.splitlines() if re.match(r'\| K[1-7] \|', line)]
metrics_preserved = len(indicator_rows(new_baseline)) == 7 and indicator_rows(old_baseline) == indicator_rows(new_baseline)
if not metrics_preserved:
    errors.append('Original K1-K7 changed')

docs = [ROOT / 'PROJECT_STATUS.md', baseline, ROOT / 'design/DESIGN_BRIEF.md', ROOT / 'README.md', ROOT / 'RUNNING.md', ROOT / '02_人体骨骼晶格肌肉方案/README.md', D / 'SHOULDER_DESIGN.md', D / 's1_GH_assembly_review/compact_recessed_keepers/README.md', Y / 'README.md', Y / 'print_inspection/README.md', Y / 'vent_experiment/README.md', ROOT / 'design/TPU_PRINTING_REFERENCE.md', ROOT / 'design/PHYSICAL_PROTOTYPE_STAGE.md', ROOT / 'docs/PROJECT_OVERVIEW.md']
for p in [ACTIVE/'README.md', Y/'transition_stock_comparison/README.md', Y/'vent_experiment/outward_start_revision/README.md', ACTIVE/'free_terminal_kinematics/README.md', ACTIVE/'free_strand_routes/README.md', ACTIVE/'volume_mesh/README.md', D/'SHOULDER_DESIGN_v0.9.16_ARCHIVE.md', ROOT/'reviews/S1_2026-11-03/document_optimization/README.md', ROOT/'reviews/S1_2026-11-03/indexed_boundary_meshing/README.md', D/'s1_tpu_channel_rebuild_revision/infraspinatus/radius_0.9/README.md']:
    if not p.exists():continue
    if p not in docs:docs.append(p)
link_count = 0
for doc in docs:
    text = doc.read_text(encoding='utf-8')
    # Historical brief kept verbatim; validate the active entry only.
    if doc.name == 'DESIGN_BRIEF.md':
        text = text.split('## 历史记录阅读说明', 1)[0]
    for target in re.findall(r'!?\[[^\]]*\]\(([^\n)]+)\)', text):
        target = unquote(target.strip('<>').split('#', 1)[0])
        if not target or '://' in target or target.startswith('mailto:'):
            continue
        link_count += 1
        if not (doc.parent / target).exists():
            errors.append(f'Broken local link: {doc.relative_to(ROOT)} -> {target}')

checks = {
    'scope': 'Document/source consistency, with existing geometric reports inspected; not a fresh geometric recomputation or a physical test.',
    'selected_source_count': len(report['sources']),
    'selected_source_hashes_match': not any('Source hash:' in e for e in errors),
    'changed_source_keys': changed,
    'unchanged_source_count': len(report['sources']) - len(changed),
    'local_links_checked': link_count,
    'historical_brief_bytes_preserved': history_preserved,
    'original_K1_K7_rows_preserved': metrics_preserved,
    'existing_report_evidence': {
        'zero_pose_pairs': report['pair_count_checked'],
        'paper_Y_frames_reported': report.get('paper_Y_frame_checks_pass'),
        'unchanged_lattice_graph_reference': report.get('frame_and_lattice_graph_reference'),
        'all_axial_slabs_supported': report['all_axial_slabs_have_supported_cells'],
        'no_rigid_overlap_over_threshold': report['no_rigid_overlap_over_threshold'],
        'no_wire_overlap_over_threshold': report['no_wire_overlap_over_threshold'],
        'stage_pass': report['stage_pass'],
        'actual_drive_range_deg': report['actual_drive_range_deg'],
    },
    'errors': errors,
    'document_sha256': {str(p.relative_to(ROOT)): sha(p) for p in docs},
}
(HERE / 'verification_checks.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in checks.items() if k != 'document_sha256'}, ensure_ascii=True, indent=2))
if errors:
    raise SystemExit(1)
