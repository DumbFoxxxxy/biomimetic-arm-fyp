"""One-time source migration; preserves model meshes and evidence JSON bytes."""
import ast
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BOOT='''
# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
'''


def run():
    changed=[]
    for folder in ['mujoco_arm','knowledge_base']:
        for p in (ROOT/folder).rglob('*.py'):
            if '__pycache__' in p.parts:continue
            original=p.read_text(encoding='utf-8-sig');text=original
            text=re.sub(r'from pathlib import Path\b(?! as)', 'from project_runtime import WorkspacePath as Path',text)
            text=text.replace('mujoco.MjModel.from_xml_path(', '_load_project_model(')
            text=text.replace('mujoco.MjModel.from_xml_string(', '_load_project_xml(')
            if text==original:continue
            if '_project_root = next(' not in text:
                tree=ast.parse(text);position=0
                for node in tree.body:
                    if isinstance(node,ast.Expr) and isinstance(node.value,ast.Constant) and isinstance(node.value.value,str) and position==0:position=node.end_lineno
                    elif isinstance(node,ast.ImportFrom) and node.module=='__future__':position=node.end_lineno
                    else:break
                lines=text.splitlines(keepends=True)
                imports='from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml\n'
                text=''.join(lines[:position])+BOOT+imports+''.join(lines[position:])
            ast.parse(text,filename=str(p))
            p.write_text(text,encoding='utf8');changed.append(str(p.relative_to(ROOT)))
    (ROOT/'整理记录/路径改造/changed_sources.json').write_text(json.dumps(changed,ensure_ascii=False,indent=2),encoding='utf8')
    print('Updated source modules:',len(changed))


if __name__=='__main__':run()
