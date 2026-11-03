"""Reuse unchanged module meshes only after exact per-module route comparison.

Rebinds metadata, not geometry or validation scope. The new whole assembly
must still be checked, particularly its changed neighbouring interfaces.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--routes',type=Path,required=True)
    ap.add_argument('--report',type=Path,action='append',required=True)
    ap.add_argument('--tag',default='upper_master');args=ap.parse_args()
    if not args.tag.replace('_','').isalnum():raise ValueError('Invalid output tag')
    new=json.loads(args.routes.read_text(encoding='utf8'));digest=hashlib.sha256(args.routes.read_bytes()).hexdigest()
    for source in args.report:
        original=json.loads(source.read_text(encoding='utf8'));module=original['module']
        old_routes=D.frames.PROJECT_ROOT/original['master_routes']
        if hashlib.sha256(old_routes.read_bytes()).hexdigest()!=original['master_sha256']:raise ValueError('Old routes changed')
        old=json.loads(old_routes.read_text(encoding='utf8'))
        old_rows=[r for r in old['routes'] if r['module']==module];new_rows=[r for r in new['routes'] if r['module']==module]
        if not old_rows or old_rows!=new_rows:raise ValueError('Own routes changed; regenerate '+module)
        for row in original['parts'].values():
            if hashlib.sha256((D.frames.PROJECT_ROOT/row['file']).read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError('Geometry changed '+module)
        result=copy.deepcopy(original)
        result.update(master_routes=str(args.routes.relative_to(D.frames.PROJECT_ROOT)),master_sha256=digest,
            own_route_records_identical=True,route_rebinding_parent=str(source.relative_to(D.frames.PROJECT_ROOT)),
            route_rebinding_parent_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),stage_pass=False)
        for key in ['unchanged_neighbour_intersections_mm3','repaired_TPU_bone_and_neighbour_intersections_mm3']:
            if key in result:result['parent_'+key]=result.pop(key)
        folder=D.DEST/('unified_'+module+'_'+args.tag);folder.mkdir(exist_ok=True)
        (folder/'check.json').write_text(json.dumps(result,indent=2),encoding='utf8')
        print('Exact unchanged route/mesh reuse',module,flush=True)


if __name__=='__main__':main()
