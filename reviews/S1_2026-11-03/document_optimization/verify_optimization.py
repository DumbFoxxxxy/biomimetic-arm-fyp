"""Verify documentation edits without regenerating geometry or mechanics."""
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
D = ROOT / '02_人体骨骼晶格肌肉方案/输出/shoulder_design'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))

def text(path):
    return path.read_text(encoding='utf-8')

errors = []
before = read(HERE / 'before_SHA256.json')
for rel, expected in before.items():
    if sha(HERE / 'before' / rel) != expected:
        errors.append('Snapshot changed: ' + rel)

selection_path = ROOT / 'design/S1_CURRENT_CANDIDATE.json'
selection = read(selection_path)
pointer_unchanged = sha(selection_path) == before['design/S1_CURRENT_CANDIDATE.json']
if not pointer_unchanged:
    errors.append('Candidate selection changed')
report_path = ROOT / selection['assembly_report']
if sha(report_path) != selection['assembly_sha256']:
    errors.append('Selected assembly changed')
report = read(report_path)
for key, source in report['sources'].items():
    if sha(ROOT / source['file']) != source['sha256']:
        errors.append('Selected source changed: ' + key)
if selection['stage_pass'] or report['stage_pass'] or selection['actual_drive_range_deg'] is not None:
    errors.append('Unexpected performance claim')

old_path = HERE / 'before' / D.relative_to(ROOT) / 'SHOULDER_DESIGN.md'
archive = D / 'SHOULDER_DESIGN_v0.9.16_ARCHIVE.md'
archive_exact = archive.read_bytes() == old_path.read_bytes()
if not archive_exact:
    errors.append('Technical archive differs from original bytes')
new = text(D / 'SHOULDER_DESIGN.md')
old = text(old_path)
rows = lambda s: [line for line in s.splitlines() if line.startswith('| K')]
metrics_exact = rows(old) == rows(new) and len(rows(new)) == 7
if not metrics_exact:
    errors.append('Technical K1-K7 rows changed')

parameter_rows = []
for line in old.splitlines():
    if any(line.startswith('| ' + label + ' |') for label in [
        '原近端盖内面至骨前面，两轴', '两轴处原骨片厚度', 'TPU夹持端垫厚度',
        '穿骨螺钉候选', '孔心距／穿孔直径', '名义X/Y/Z尺寸',
        '单元／轴向分层／表面面片', '实际轴向站距', '实际横向X投影间距范围',
        '实际横向Z投影间距范围', 'TPU体积', '冈下肌肌肉侧', '冈下肌骨侧',
        '肩胛下肌肌肉侧', '肩胛下肌骨侧',
    ]):
        parameter_rows.append(line)
        expected = line.replace('| TPU体积 |', '| Y1父级TPU体积 |')
        if expected not in new:
            errors.append('Parameter row missing: ' + line)

marker = '## 历史记录阅读说明'.encode('utf-8')
brief_before = (HERE / 'before/design/DESIGN_BRIEF.md').read_bytes()
brief_now = (ROOT / 'design/DESIGN_BRIEF.md').read_bytes()
history_exact = brief_before.split(marker, 1)[1] == brief_now.split(marker, 1)[1]
if not history_exact:
    errors.append('Historical brief changed')
baseline_metrics_exact = rows(text(HERE / 'before/design/CURRENT_BASELINE.md')) == rows(text(ROOT / 'design/CURRENT_BASELINE.md'))
if not baseline_metrics_exact:
    errors.append('Baseline K1-K7 rows changed')

branches = {}
for rel in ['s1_tpu_distal_surface_revision/junction_repair_trials.json', 's1_tpu_surface_global_revision/global_repair_trials.json']:
    p = D / rel
    data = read(p)
    if data['adopted_as_current'] or data['stage_pass'] or data['actual_drive_range_deg'] is not None:
        errors.append('Experimental branch status changed: ' + rel)
    if data['parent_sha256'] != sha(report_path):
        errors.append('Experimental parent mismatch: ' + rel)
    branch = {'report_sha256': sha(p), 'adopted_as_current': False, 'stage_pass': False}
    if 'results' in data:
        exports = []
        for module in data['results']:
            for e in module['exports']:
                mesh = e['mesh']
                if sha(ROOT / mesh['file']) != mesh['sha256']:
                    errors.append('Experimental export hash: ' + mesh['file'])
                count = e['surface_audit']['count']
                if count <= 0:
                    errors.append('Unexpected clean export; review prose: ' + mesh['file'])
                exports.append({'module': module['module'], 'format': e['kind'], 'recorded_detected_crossings': count, 'mesh_sha256': mesh['sha256']})
        branch['recorded_export_diagnostics_not_recomputed'] = exports
    branches[rel] = branch

result = {
    'scope': 'Documentation/source consistency only; no new geometry, meshing, equilibrium or physical validation.',
    'selected_source_count': len(report['sources']),
    'candidate_pointer_byte_identical': pointer_unchanged,
    'selected_assembly_sha256': sha(report_path),
    'technical_archive_byte_identical': archive_exact,
    'technical_lines_before': len(old.splitlines()),
    'technical_lines_after': len(new.splitlines()),
    'technical_K1_K7_rows_byte_content_preserved': metrics_exact,
    'baseline_K1_K7_preserved': baseline_metrics_exact,
    'historical_task_brief_bytes_preserved': history_exact,
    'original_parameter_rows_checked': len(parameter_rows),
    'experimental_report_status_only': branches,
    'stage_pass': False,
    'actual_drive_range_deg': None,
    'errors': errors,
}
canonical_path = ROOT / 'reviews/S1_2026-11-02/Y_muscle_pair/verification_checks.json'
canonical = read(canonical_path)
if canonical['errors']:
    errors.append('Canonical document/source checks failed')
result['local_links_checked_by_canonical_verifier'] = canonical['local_links_checked']
result['document_sha256'] = canonical['document_sha256']
(HERE / 'verification_checks.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in result.items() if k not in ['experimental_report_status_only', 'document_sha256']}, ensure_ascii=True, indent=2))
if errors:
    raise SystemExit(1)
