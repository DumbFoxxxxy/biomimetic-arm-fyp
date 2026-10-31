"""Read-only model/source audit plus a small, real alternate-root viewer copy."""
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from project_runtime import project_path,load_mujoco_model


def run():
    report={'scenes':[],'source_hashes':[],'parse_failures':[]}
    for p in (ROOT/'mujoco_arm').rglob('*.py'):
        try:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
        except Exception as e:report['parse_failures'].append({'file':str(p),'error':str(e)})
    for rel in ['mujoco_arm/arm_tendon_demo.xml',
        'mujoco_arm/exports/skeleton_from_v2/V2-S1/skeleton.xml',
        'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/continuous_channel_candidate/neutral_view.xml',
        'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/rigid_hinges.xml']:
        m=load_mujoco_model(rel);report['scenes'].append({'file':str(project_path(rel).relative_to(ROOT)),'geoms':m.ngeom,'nq':m.nq})
    out=project_path('mujoco_arm/exports/anatomical_concept')
    for rel,digest in json.loads((out/'elbow_handoff_provenance.json').read_text(encoding='utf8'))['source_sha256'].items():
        p=project_path(rel);actual=hashlib.sha256(p.read_bytes()).hexdigest()
        report['source_hashes'].append({'file':str(p.relative_to(ROOT)),'unchanged':actual==digest})
    # A fresh ordinary directory: no symlinks/junctions, no copy of original root.
    clone=ROOT/'tmp/路径验证 克隆目录'
    clone.mkdir(parents=True,exist_ok=True)
    shutil.copytree(ROOT/'project_runtime',clone/'project_runtime',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__'))
    script='mujoco_arm/anatomical_concept/view_continuous_elbow.py'
    (clone/script).parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/script,clone/script)
    base=out/'direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
    xmls=[base/'neutral_view.xml',base/'continuous_channel_candidate/neutral_view.xml']
    for xml in xmls:
        dest=clone/xml.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(xml,dest)
        for mesh in ET.fromstring(xml.read_text(encoding='utf-8-sig')).findall('.//asset/mesh'):
            src=project_path(mesh.attrib['file']);dst=clone/src.relative_to(ROOT)
            dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    p=subprocess.run([sys.executable,str(clone/script),'--check'],cwd=ROOT/'tools',stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf8',errors='replace')
    report['alternate_unicode_clone']={'root':str(clone),'exit_code':p.returncode,'output':p.stdout.strip(),'uses_junctions':False}
    assert p.returncode==0,p.stdout
    assert not report['parse_failures'] and all(x['unchanged'] for x in report['source_hashes'])
    (ROOT/'整理记录/路径改造/verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'scene_counts':[(s['geoms'],s['nq']) for s in report['scenes']],'source_hashes_unchanged':True,'alternate_unicode_clone_passed':True,'python_AST_passed':True}))


if __name__=='__main__':run()
