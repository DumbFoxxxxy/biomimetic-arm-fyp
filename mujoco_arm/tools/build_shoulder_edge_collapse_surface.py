"""Explicit local short-edge collapses for complete TPU analysis candidates.

Retains a connected closed indexed surface using the edge link condition,
checks surviving-face orientation, and bounds each merged vertex cluster's
displacement. This is a reported geometry approximation, not a manufacturing
or mechanical release. Original source files and selection remain unchanged.
"""
import argparse
import heapq
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom


def collapse_short_edges(vertices, faces, tolerance):
    V=vertices.copy(); F=faces.copy(); active=np.ones(len(F),dtype=bool)
    incident=[set() for _ in V]
    for index,face in enumerate(F):
        for vertex in face:incident[vertex].add(index)
    clusters={i:[i] for i in range(len(V))}
    edges=np.unique(np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0)
    queue=[]
    for a,b in edges:
        length=float(np.linalg.norm(V[a]-V[b]))
        if length<=tolerance:heapq.heappush(queue,(length,int(a),int(b)))
    operations=[]; rejected={'link':0,'orientation':0,'cluster_bound':0}

    def neighbours(vertex):
        return set(F[list(incident[vertex])].ravel())-{vertex}

    while queue:
        _,a,b=heapq.heappop(queue)
        if not incident[a] or not incident[b] or np.linalg.norm(V[a]-V[b])>tolerance:continue
        shared=incident[a]&incident[b]
        if len(shared)!=2:continue
        opposites=set(F[list(shared)].ravel())-{a,b}
        if neighbours(a)&neighbours(b)!=opposites:
            rejected['link']+=1;continue
        members=clusters[a]+clusters[b]; position=(V[a]+V[b])/2
        max_move=float(np.linalg.norm(vertices[members]-position,axis=1).max())
        if max_move>tolerance:
            rejected['cluster_bound']+=1;continue
        affected=incident[a]|incident[b]; surviving=sorted(affected-shared)
        before=V[F[surviving]]; proposed=F[surviving].copy();proposed[proposed==b]=a
        after=V[proposed].copy();after[proposed==a]=position
        n0=np.cross(before[:,1]-before[:,0],before[:,2]-before[:,0])
        n1=np.cross(after[:,1]-after[:,0],after[:,2]-after[:,0])
        guard=np.linalg.norm(n0,axis=1)>1e-12
        if np.any(np.linalg.norm(n1,axis=1)<=1e-14) or np.any(np.einsum('ij,ij->i',n0[guard],n1[guard])<=0):
            rejected['orientation']+=1;continue
        for index in affected:
            for vertex in F[index]:incident[vertex].discard(index)
        for index in shared:active[index]=False
        V[a]=position;F[surviving]=proposed;clusters[a]=members;clusters.pop(b)
        for index in surviving:
            for vertex in F[index]:incident[vertex].add(index)
        operations.append(dict(retained_vertex=a,removed_vertex=b,new_position_source_mm=position.tolist(),
            original_cluster_vertex_count=len(members),maximum_original_cluster_displacement_mm=max_move))
        for vertex in neighbours(a):
            low,high=sorted([a,int(vertex)]);length=float(np.linalg.norm(V[low]-V[high]))
            if length<=tolerance:heapq.heappush(queue,(length,low,high))
    used=np.unique(F[active]); mapping=np.full(len(V),-1,dtype=np.int64);mapping[used]=np.arange(len(used))
    return V[used],mapping[F[active]],dict(tolerance_mm=tolerance,collapses=len(operations),operations=operations,
        rejected_candidates=rejected,maximum_cluster_displacement_mm=max((r['maximum_original_cluster_displacement_mm'] for r in operations),default=0.))


def run(module,tolerance):
    if not 0<tolerance<=.02:raise ValueError('Explicit research approximation must be 0 < tolerance <= 0.02mm')
    selection=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
    parent_file=A.ROOT/selection['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=selection['assembly_sha256']:raise ValueError('Stale selected parent')
    source=parent['sources']['tpu:'+module];original=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])[0]
    V,F=mesh_geom.manifold_mesh_to_arrays(original.to_mesh64());V,F,history=collapse_short_edges(V,F,tolerance)
    _,boundary,nonmanifold=mesh_geom.edge_census(F)
    if boundary or nonmanifold:raise ValueError('Short-edge collapse failed topology condition')
    actual=A.I._manifold().Manifold(A.I._manifold().Mesh64(np.ascontiguousarray(V),np.ascontiguousarray(F,dtype=np.uint64)))
    if str(actual.status()).split('.')[-1]!='NoError' or actual.volume()<=0:raise ValueError('Collapsed surface is invalid')
    dest=A.DESIGN/'s1_tpu_edge_collapse_revision'/module/('epsilon_'+f'{tolerance:g}');dest.mkdir(parents=True,exist_ok=True)
    path=dest/'TPU_analysis_surface.obj'
    with path.open('w',encoding='utf8',newline='\n') as stream:
        stream.write('# Complete TPU research analysis approximation; source coordinates mm; no weld\n')
        for v in V:stream.write('v %.17g %.17g %.17g\n'%tuple(v))
        for f in F+1:stream.write('f %d %d %d\n'%tuple(f))
    vv,ff=mesh_geom.load_obj(str(path))
    if not np.array_equal(V,vv) or not np.array_equal(F,ff):raise ValueError('Indexed OBJ roundtrip changed coordinates/indices')
    tri=V[F];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    row=dict(file=path.relative_to(A.ROOT).as_posix(),sha256=A.hashed(path),translation_source_mm=[0.,0.,0.],
        import_convention='indexed_float64_obj_no_bucket_weld',boundary_edges=boundary,nonmanifold_edges=nonmanifold,
        volume_mm3=float(actual.volume()))
    sources=dict(parent['sources']);sources['tpu:'+module]=row
    report=dict(scope=__doc__,module=module,source=source,output=row,script_sha256=A.hashed(Path(__file__).resolve()),
        parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),collapse_history=history,
        triangles=len(F),minimum_triangle_area_mm2=float(areas.min()),positive_components=sum(p.volume()>1e-6 for p in actual.decompose()),
        added_material_signed_mm3=float((actual-original).volume()),removed_material_signed_mm3=float((original-actual).volume()),
        source_bytes_unchanged=A.hashed(A.ROOT/source['file'])==source['sha256'],
        geometric_self_intersection_free_verified=False,interfaces_and_channels_reverified=False,
        manufacturing_approximation_adopted=False,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (dest/'surface_revision.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    assembly=dict(scope=__doc__,sources=sources,parent_report=report['parent_report'],parent_sha256=report['parent_sha256'],
        source_revision_report=(dest/'surface_revision.json').relative_to(A.ROOT).as_posix(),stage_pass=False,actual_drive_range_deg=None)
    (dest/'candidate_assembly.json').write_text(json.dumps(assembly,indent=2)+'\n',encoding='utf8')
    print(module,tolerance,'collapses',history['collapses'],'faces',len(F),'min area',areas.min(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--module',required=True,choices=['infraspinatus','subscapularis'])
    parser.add_argument('--tolerance',type=float,default=.001);args=parser.parse_args();run(args.module,args.tolerance)
