"""Export exact source poses for the verified conditional bundle corridor.

No interpolation of muscle shape or claimed physical actuation. The explicit
later-stage inventory is absent only from this staged view, never the design.
"""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf8'))
def save(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')


def main():
    source=DESIGN/'s1_deltoid_support_draft/assembly_check.json';r=read(source)
    proof_file=source.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(source) and proof['conditional_bundle_placement_geometry_pass']
    assert proof['verifier_sha256']==sha(ROOT/'tools/verify_s1_deltoid_support_draft.py')
    for row in r['sources'].values():assert sha(ROOT/row['file'])==row['sha256']
    destination=source.parent/'entry_demo';destination.mkdir(exist_ok=True)
    present=set(proof['fixed_keys']+proof['moving_keys']);direction=proof['direction_source']
    for offset in [0,80]:
        sources={k:dict(row) for k,row in r['sources'].items() if k in present}
        for key in proof['moving_keys']:
            old=sources[key]['translation_source_mm'] or [0.,0.,0.]
            sources[key]['translation_source_mm']=[old[i]+direction[i]*offset for i in range(3)]
        save(destination/f'pose_{offset}mm.json',dict(scope=__doc__,sources=sources,n_instances=len(sources),
            source_report=source.relative_to(ROOT).as_posix(),source_report_sha256=sha(source),verification_sha256=sha(proof_file),
            script_sha256=sha(Path(__file__)),assembly_offset_mm=offset,moving_keys=proof['moving_keys'],
            omitted_later_keys=proof['later_installation_keys'],temporary_fixture_designed=False,
            whole_assembly_verified=False,driven_motion_verified=False,stage_pass=False))
        print('pose',offset,'present',len(sources),'pending',len(proof['later_installation_keys']),flush=True)


if __name__=='__main__':main()
