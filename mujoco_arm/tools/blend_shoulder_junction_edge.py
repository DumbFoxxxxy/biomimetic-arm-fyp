"""Add a measured local TPU blend at a self-touching lattice junction edge.

Adds actual material, keeps strict topology checks, then rechecks complete
steel routes and installed caps. No acceptance threshold is relaxed.
"""
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import mesh_geom
import interference as I
import check_shoulder_assembly as C
import build_shoulder_tendon_segments as T
from build_shoulder_unified_muscle import local_export


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path,required=True)
    ap.add_argument('--radius',type=float,default=.2);args=ap.parse_args()
    if not .1<=args.radius<=1.:raise ValueError('Inspect before using a blend larger than the 1mm lattice radius')
    meta=json.loads(args.report.read_text(encoding='utf8'));row=meta['parts']['TPU']
    path=D.frames.PROJECT_ROOT/row['file']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed TPU source')
    V,F=mesh_geom.read_stl(str(path));W,G=mesh_geom.weld(V,F)
    edges=np.sort(np.vstack([G[:,[0,1]],G[:,[1,2]],G[:,[2,0]]]),axis=1);edge,count=np.unique(edges,axis=0,return_counts=True)
    bad=edge[count>2]
    if len(bad)!=1:raise ValueError('Expected one measured junction edge; inspect changed topology')
    origin=np.array(row['translation_source_mm']);points=W[bad[0]]+origin
    source=D.frames.PROJECT_ROOT/meta['master_routes']
    if hashlib.sha256(source.read_bytes()).hexdigest()!=meta['master_sha256']:raise ValueError('Changed master routes')
    routes=[r for r in json.loads(source.read_text(encoding='utf8'))['routes'] if r['module']==meta['module']]
    original=C.load_stl(path).translate(origin)
    blend=D.outer_envelope(points,args.radius)
    bones={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()}
    def restore(r):return C.load_stl(D.frames.PROJECT_ROOT/r['file']).translate(r['translation_source_mm'])
    caps={k:restore(r) for k,r in meta['parts'].items() if k!='TPU'}
    for key,s in {**bones,**caps}.items():
        if I.intersection_volume(blend,s)>.001:raise ValueError('Blend touches installed part '+key)
    lumens=I._manifold().Manifold.batch_boolean([D.outer_envelope(np.array(r['full_wire_path_source_mm']),.86) for r in routes],I._manifold().OpType.Add)
    candidate=(original+blend)-lumens
    tag=str(args.radius).replace('.','p')
    folder=D.DEST/('unified_'+meta['module']+'_upper_blend_r'+tag);folder.mkdir(exist_ok=True)
    solid,part=local_export(candidate,folder/'TPU_muscle.stl',origin)
    result=copy.deepcopy(meta);result['parts']['TPU']=part
    result.update(topology_blend=dict(radius_mm=args.radius,edge_source_mm=points.tolist(),
        added_material_mm3=float((candidate-original).volume()),removed_material_mm3=float((original-candidate).volume()),
        parent_report=str(args.report.relative_to(D.frames.PROJECT_ROOT)),
        parent_sha256=hashlib.sha256(args.report.read_bytes()).hexdigest()),stage_pass=False)
    for key in ['free_TPU_wall_guard_outside_assigned_cap_housings','missing_0_75mm_wall_guard_mm3','unchanged_neighbour_intersections_mm3']:
        if key in result:result['parent_'+key]=result.pop(key)
    result['wall_guard_after_blend_verified']=False
    parts={'TPU':solid,**caps}
    result['wire_material_intersections_mm3']={r['head']:{k:I.intersection_volume(
        D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3),s) for k,s in parts.items()} for r in routes}
    result['material_pair_intersections_mm3']=[dict(a=a,b=b,intersection_mm3=I.intersection_volume(parts[a],parts[b])) for a,b in itertools.combinations(parts,2)]
    result['all_part_bone_intersections_mm3']={k:{b:I.intersection_volume(s,v) for b,v in bones.items()} for k,s in parts.items()}
    result['bone_intersections_mm3']=result['all_part_bone_intersections_mm3']['TPU']
    (folder/'check.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print('Local blend mesh',part['mesh_checks_pass'],'edges',part['boundary_edges'],part['nonmanifold_edges'],
        'material added',result['topology_blend']['added_material_mm3'],flush=True)


if __name__=='__main__':main()
