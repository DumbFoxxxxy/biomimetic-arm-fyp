"""Bounded local Laplacian fairing of degenerate TPU analysis facets.

Changes only a recorded neighbourhood of poor triangles in a short-edge
candidate. Maximum vertex displacement is explicit; this is a geometric
approximation requiring interface/channel checks, not a printable release.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom


def triangle_stats(V,F):
    t=V[F];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);area=np.linalg.norm(normal,axis=1)/2
    edge2=sum(np.sum((t[:,i]-t[:,j])**2,axis=1) for i,j in [(0,1),(1,2),(2,0)])
    return area,4*np.sqrt(3)*area/np.maximum(edge2,1e-300),normal


def run(module,maximum_move):
    if not 0<maximum_move<=.02:raise ValueError('Research fairing bounded to at most 0.02mm per vertex')
    selection=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
    parent_file=A.ROOT/selection['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=selection['assembly_sha256']:raise ValueError('Stale selected parent')
    input_folder=A.DESIGN/'s1_tpu_edge_collapse_revision'/module/'epsilon_0.001'
    input_report=input_folder/'surface_revision.json';r=json.loads(input_report.read_text(encoding='utf8'))
    mesh=r['output'];path=A.ROOT/mesh['file']
    if A.hashed(path)!=mesh['sha256'] or r['parent_sha256']!=A.hashed(parent_file):raise ValueError('Stale fairing input')
    initial,F=mesh_geom.load_obj(str(path));V=initial.copy();area,quality,n0=triangle_stats(V,F)
    bad=(quality<1e-4)|(area<1e-6)
    active=np.zeros(len(V),dtype=bool);active[np.unique(F[bad])]=True
    edges=np.unique(np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0)
    ring=edges[np.any(active[edges],axis=1)];active[np.unique(ring)]=True
    degree=np.bincount(edges.ravel(),minlength=len(V));history=[]
    for iteration in range(12):
        total=np.zeros_like(V);np.add.at(total,edges[:,0],V[edges[:,1]]);np.add.at(total,edges[:,1],V[edges[:,0]])
        proposed=V.copy();proposed[active]+=.25*(total[active]/degree[active,None]-V[active])
        delta=proposed-initial;distance=np.linalg.norm(delta,axis=1)
        scale=np.minimum(1.,maximum_move/np.maximum(distance,1e-300));V=initial+delta*scale[:,None]
        area,q,n1=triangle_stats(V,F)
        history.append(dict(iteration=iteration+1,triangles_below_quality_1e_4=int(np.sum(q<1e-4)),
            triangles_below_area_1e_6_mm2=int(np.sum(area<1e-6)),minimum_area_mm2=float(area.min()),
            maximum_vertex_displacement_mm=float(np.linalg.norm(V-initial,axis=1).max())))
    # Preserve precisely the exported vertex positions; native CSG is diagnostic.
    solid=A.I._manifold().Manifold(A.I._manifold().Mesh64(np.ascontiguousarray(V),np.ascontiguousarray(F,dtype=np.uint64)))
    if solid.volume()<=0 or str(solid.status()).split('.')[-1]!='NoError':raise ValueError('Fairing broke indexed topology')
    original_row=parent['sources']['tpu:'+module]
    original=A.load_checked(A.ROOT/original_row['file'],original_row['translation_source_mm'],original_row['sha256'])[0]
    dest=A.DESIGN/'s1_tpu_local_fairing_revision'/module/('move_'+f'{maximum_move:g}');dest.mkdir(parents=True,exist_ok=True)
    path=dest/'TPU_analysis_surface.obj'
    with path.open('w',encoding='utf8',newline='\n') as stream:
        stream.write('# Bounded local fairing research approximation; source mm; no welding\n')
        for v in V:stream.write('v %.17g %.17g %.17g\n'%tuple(v))
        for f in F+1:stream.write('f %d %d %d\n'%tuple(f))
    vv,ff=mesh_geom.load_obj(str(path))
    if not np.array_equal(V,vv) or not np.array_equal(F,ff):raise ValueError('Float64 OBJ roundtrip changed')
    _,boundary,nm=mesh_geom.edge_census(F)
    output=dict(file=path.relative_to(A.ROOT).as_posix(),sha256=A.hashed(path),translation_source_mm=[0.,0.,0.],
        import_convention='indexed_float64_obj_no_bucket_weld',boundary_edges=boundary,nonmanifold_edges=nm,volume_mm3=float(solid.volume()))
    report=dict(scope=__doc__,module=module,script_sha256=A.hashed(Path(__file__).resolve()),
        parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
        input_revision_report=input_report.relative_to(A.ROOT).as_posix(),input_revision_sha256=A.hashed(input_report),
        input_mesh=mesh,output=output,active_vertices=int(active.sum()),maximum_allowed_vertex_displacement_mm=maximum_move,
        iterations=history,faces_with_reversed_normal_from_input=int(np.sum(np.einsum('ij,ij->i',n0,n1)<0)),
        added_material_signed_mm3=float((solid-original).volume()),removed_material_signed_mm3=float((original-solid).volume()),
        source_bytes_unchanged=A.hashed(A.ROOT/original_row['file'])==original_row['sha256'],
        geometric_self_intersection_free_verified=False,interfaces_and_channels_reverified=False,
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (dest/'surface_revision.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    sources=dict(parent['sources']);sources['tpu:'+module]=output
    assembly=dict(scope=__doc__,sources=sources,parent_report=report['parent_report'],parent_sha256=report['parent_sha256'],
        source_revision_report=(dest/'surface_revision.json').relative_to(A.ROOT).as_posix(),stage_pass=False,actual_drive_range_deg=None)
    (dest/'candidate_assembly.json').write_text(json.dumps(assembly,indent=2)+'\n',encoding='utf8')
    print(module,'active',active.sum(),'bad after',history[-1],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--module',choices=['infraspinatus','subscapularis'],required=True)
    p.add_argument('--maximum-move',type=float,default=.02);args=p.parse_args();run(args.module,args.maximum_move)
