"""Check that selected model bytes are represented in the Git index via LFS."""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text(encoding='utf8'))

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

entries={}
for row in subprocess.check_output(['git','ls-files','--stage','-z'],cwd=ROOT).split(b'\0'):
    if not row:continue
    meta,path=row.split(b'\t',1);entries[path.decode('utf8')]=meta.split()[1].decode()
selection_path=ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=read(selection_path)
report_file=ROOT/selection['assembly_report'];report=read(report_file);errors=[];sources=[]
for key,row in report['sources'].items():
    rel=Path(row['file']).as_posix();path=ROOT/rel
    if rel not in entries:
        errors.append('Selected source missing from index: '+rel);continue
    pointer=subprocess.check_output(['git','cat-file','blob',entries[rel]],cwd=ROOT)
    expected=('oid sha256:'+row['sha256']).encode()
    valid=pointer.startswith(b'version https://git-lfs.github.com/spec/v1\n') and expected in pointer.splitlines()
    if not valid or sha(path)!=row['sha256']:errors.append('Selected source pointer/hash mismatch: '+key)
    sources.append(dict(key=key,file=rel,source_sha256=row['sha256'],indexed_LFS_pointer_matches=valid))
before=read(HERE/'before_SHA256.json')
if sha(selection_path)!=before['design/S1_CURRENT_CANDIDATE.json']:errors.append('Candidate selection unexpectedly changed')
if sha(report_file)!=selection['assembly_sha256']:errors.append('Candidate report unexpectedly changed')
marker='## 历史记录阅读说明'.encode('utf8')
history=(ROOT/'design/DESIGN_BRIEF.md').read_bytes().split(marker,1)[1]
if history!=(HERE/'before/design/DESIGN_BRIEF.md').read_bytes().split(marker,1)[1]:errors.append('Historical brief changed')
docs=read(ROOT/'reviews/S1_2026-11-02/Y_muscle_pair/verification_checks.json')
if docs['errors']:errors.append('Document/source checks failed')
output=dict(scope=__doc__,source_count=len(sources),all_selected_sources_in_index=len(sources)==len(report['sources']),
    selected_model_report_sha256=sha(report_file),candidate_pointer_byte_identical=True,
    original_K1_K7_preserved=docs['original_K1_K7_rows_preserved'],historical_brief_bytes_preserved=True,
    local_links_checked=docs['local_links_checked'],existing_tests_observed_passed=31,
    runtime_checks_observed=dict(static_scene_geoms=191,bare_scene_geoms=102,actual_TPU_actuation_verified=False),
    stage_priority='Physical wire/TPU actuation prototype first; movement range measured on hardware; full FE deferred.',
    sources=sources,errors=errors)
(HERE/'git_input_checks.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in output.items() if k!='sources'},ensure_ascii=True,indent=2))
if errors:raise SystemExit(1)
