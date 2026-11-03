"""Explicit bore geometry revision, independent of the selected S1-Y2 model.

Re-cuts a slightly larger channel along the unchanged full wire path with
positively overlapping sphere/cylinder joints. Intentional bore removal is
reported separately from unintended loss; no collision tolerance is relaxed.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
import build_shoulder_tpu as T
from build_shoulder_shaft_collar import rod


def run(module,radius):
    if not .86<radius<=.95:raise ValueError('Explicit bore research radius must be >0.86 and <=0.95mm')
    choice=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
    parent_file=A.ROOT/choice['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=choice['assembly_sha256']:raise ValueError('Stale selected input')
    row=parent['sources']['tpu:'+module];old=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])[0]
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    route=next(r for r in json.loads(route_file.read_text(encoding='utf8'))['routes'] if r['module']==module)
    path=np.array(route['full_wire_path_source_mm']);node_radius=radius*1.01
    pieces=[T.rod(a,b,radius,segments=32) for a,b in zip(path,path[1:]) if np.linalg.norm(b-a)>1e-6]
    pieces+=[A.I._manifold().Manifold.sphere(node_radius,32).translate(p) for p in path]
    cutter=A.I._manifold().Manifold.batch_boolean(pieces,A.I._manifold().OpType.Add)
    candidate=old-cutter
    if candidate.is_empty() or str(candidate.status()).split('.')[-1]!='NoError':raise ValueError('Invalid channel revision')
    V,F=mesh_geom.manifold_mesh_to_arrays(candidate.to_mesh64());_,boundary,nm=mesh_geom.edge_census(F)
    if boundary or nm:raise ValueError('Bore revision broke indexed topology')
    folder=A.DESIGN/'s1_tpu_channel_rebuild_revision'/module/('radius_'+f'{radius:g}');folder.mkdir(parents=True,exist_ok=True)
    file=folder/'TPU_rebuilt_channel.obj'
    with file.open('w',encoding='utf8',newline='\n') as stream:
        stream.write('# Complete TPU explicit bore revision; source mm; no bucket welding\n')
        for vertex in V:stream.write('v %.17g %.17g %.17g\n'%tuple(vertex))
        for face in F+1:stream.write('f %d %d %d\n'%tuple(face))
    rv,rf=mesh_geom.load_obj(str(file))
    if not np.array_equal(V,rv) or not np.array_equal(F,rf):raise ValueError('Precision readback changed')
    guards=[]
    for port in route['ports'].values():
        c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source']);guards.append(rod(c-n*8.,c+n*3.,12.))
    end_guard=A.I._manifold().Manifold.batch_boolean(guards,A.I._manifold().OpType.Add)
    internal_paths=[route['core_path_source_mm']]+[route['connectors'][side]['path_source_mm'] for side in ['origin','distal']]
    wall_guard=A.I._manifold().Manifold.batch_boolean([A.D.outer_envelope(np.array(p),node_radius+.75) for p in internal_paths],A.I._manifold().OpType.Add)-cutter-end_guard
    expected_removed=old^cutter
    output=dict(file=file.relative_to(A.ROOT).as_posix(),sha256=A.hashed(file),translation_source_mm=[0.,0.,0.],
        import_convention='indexed_float64_obj_no_bucket_weld',boundary_edges=boundary,nonmanifold_edges=nm,volume_mm3=float(candidate.volume()))
    report=dict(scope=__doc__,module=module,source=row,output=output,script_sha256=A.hashed(Path(__file__).resolve()),
        parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),routes=route_file.relative_to(A.ROOT).as_posix(),routes_sha256=A.hashed(route_file),
        previous_nominal_channel_diameter_mm=1.72,new_cylinder_diameter_mm=radius*2,new_joint_sphere_diameter_mm=node_radius*2,
        cylinder_segments=32,sphere_segments=32,positive_joint_radius_ratio=1.01,
        wire_path_and_anchors_unchanged=True,steel_diameter_mm=.6,
        expected_explicit_bore_removal_mm3=float(expected_removed.volume()),
        removed_material_signed_mm3=float((old-candidate).volume()),added_material_signed_mm3=float((candidate-old).volume()),
        unintended_removal_outside_new_cutter_signed_mm3=float(((old-candidate)-cutter).volume()),
        missing_nominal_0_75mm_new_wall_guard_outside_end_interfaces_mm3=float((wall_guard-candidate).volume()),
        old_wall_guard_missing_mm3=float((wall_guard-old).volume()),
        positive_components=sum(p.volume()>1e-6 for p in candidate.decompose()),
        source_bytes_unchanged=A.hashed(A.ROOT/row['file'])==row['sha256'],
        interfaces_and_channels_reverified=False,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (folder/'channel_revision.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    sources=dict(parent['sources']);sources['tpu:'+module]=output
    assembly=dict(scope=__doc__,sources=sources,parent_report=report['parent_report'],parent_sha256=report['parent_sha256'],
        source_revision_report=(folder/'channel_revision.json').relative_to(A.ROOT).as_posix(),stage_pass=False,actual_drive_range_deg=None)
    (folder/'candidate_assembly.json').write_text(json.dumps(assembly,indent=2)+'\n',encoding='utf8')
    print(module,'new diameter',radius*2,'removed',report['removed_material_signed_mm3'],
        'new wall guard missing',report['missing_nominal_0_75mm_new_wall_guard_outside_end_interfaces_mm3'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--module',choices=['infraspinatus','subscapularis'],required=True)
    p.add_argument('--radius',type=float,default=.9);args=p.parse_args();run(args.module,args.radius)
