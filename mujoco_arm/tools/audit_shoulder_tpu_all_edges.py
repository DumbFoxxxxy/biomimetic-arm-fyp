"""Scan every unique boundary edge against the complete current TPU surface.

Native BVH ray hits only nominate facets. Each nonadjacent candidate is
checked independently by strict barycentric segment/facet intersection.
Floating-point tests are diagnostic, not an exact-arithmetic proof.
"""
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from diagnose_shoulder_tpu_volume_mesh import segment_facet

DEST=A.DESIGN/'s1_tpu_distal_surface_revision'

def run():
    selection=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
    parent_file=A.ROOT/selection['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=selection['assembly_sha256']:raise ValueError('Changed selected parent')
    results=[];DEST.mkdir(parents=True,exist_ok=True)
    for module in ['infraspinatus','subscapularis']:
        source=parent['sources']['tpu:'+module];solid,_,_=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])
        mesh=solid.to_mesh64();V,F=mesh_geom.manifold_mesh_to_arrays(mesh);ids=np.array(mesh.face_id,dtype=np.int64)
        lookup={}
        for i,face_id in enumerate(ids):lookup.setdefault(int(face_id),[]).append(i)
        edges=np.unique(np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0)
        started=time.monotonic();crossings={};candidate_tests=0;native_hits=0
        for edge_index,edge in enumerate(edges):
            p,q=V[edge];hits=solid.ray_cast(p,q);native_hits+=len(hits);tested=set()
            for h in hits:
                for facet in lookup.get(int(h.face_id),[]):
                    if facet in tested:continue
                    tested.add(facet)
                    if np.any(F[facet]==edge[0]) or np.any(F[facet]==edge[1]):continue
                    candidate_tests+=1;hit=segment_facet(p,q,V[F[facet]])
                    if hit is None:continue
                    point=np.array(hit['point_local_mm']);key=tuple(np.round(point,6))
                    if key not in crossings:
                        crossings[key]=dict(edge_index=edge_index,edge_vertex_indices=edge.tolist(),facet_index=int(facet),
                            edge_points_source_mm=V[edge].tolist(),facet_points_source_mm=V[F[facet]].tolist(),
                            intersection_source_mm=point.tolist(),barycentric_parameters=[hit['t'],hit['u'],hit['v']],native_face_id=int(h.face_id))
            if (edge_index+1)%20000==0:print(module,'edges',edge_index+1,'/',len(edges),'distinct intersections',len(crossings),flush=True)
        row=dict(module=module,source=source,vertices=len(V),triangles=len(F),unique_edges_scanned=len(edges),native_ray_hits=native_hits,
            independent_facet_tests=candidate_tests,distinct_intersection_count=len(crossings),intersections=list(crossings.values()),
            deduplication_coordinate_decimals=6,elapsed_s=time.monotonic()-started,
            scope='all current validated surface edges, candidate facets from native BVH; strict interiors and floating tolerances recorded in source')
        results.append(row)
        (DEST/(module+'_all_edges.json')).write_text(json.dumps(row,indent=2)+'\n',encoding='utf8')
        print(module,'finished',len(edges),'edges;',len(crossings),'distinct crossings',flush=True)
    report=dict(scope=__doc__,parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
        script_sha256=A.hashed(Path(__file__).resolve()),intersection_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/diagnose_shoulder_tpu_volume_mesh.py'),
        modules=results,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'all_edges_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
