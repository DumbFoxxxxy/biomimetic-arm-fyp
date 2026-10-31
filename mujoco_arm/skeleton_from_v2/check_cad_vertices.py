"""Independent nearest-vertex correspondence for STEP roundtrip geometry."""
import json,numpy as np
from scipy.spatial import cKDTree
from verify_joint_candidate import OUT,md
def run():
    source=json.loads((OUT/'cad_input.json').read_text());cad=json.loads((OUT/'cad_roundtrip_vertices.json').read_text());reexport=json.loads((OUT/'cad_roundtrip_meshes.json').read_text());rows=[]
    if reexport['input_cache_sha256']!=source['input_cache_sha256']:raise RuntimeError('Stale CAD reexport')
    for part in source['parts']:
        a=np.array(part['vertices_design_mm']);b=np.unique(np.array(cad[part['name']]),axis=0)
        e1=float(cKDTree(a).query(b)[0].max());e2=float(cKDTree(b).query(a)[0].max())
        inputsolid=md.Manifold(md.Mesh(np.asarray(a,dtype=np.float32),np.asarray(part['faces'],dtype=np.uint32)))
        tri=np.asarray(reexport['meshes_design_mm'][part['name']]);points,faces=np.unique(np.round(tri.reshape(-1,3),6).astype(np.float32),axis=0,return_inverse=True)
        other=md.Manifold(md.Mesh(points,np.asarray(faces.reshape(-1,3),dtype=np.uint32)))
        valid=str(other.status())=='Error.NoError'
        difference=float((inputsolid-other).volume()+(other-inputsolid).volume()) if valid else None
        limit=10 if part['name']=='thorax_geom' else 2
        rows.append({'bone':part['name'],'roundtrip_to_input_max_vertex_distance_mm':e1,'input_to_roundtrip_max_vertex_distance_mm':e2,'reexport_mesh_status':str(other.status()),'symmetric_difference_volume_mm3':difference,'difference_volume_limit_mm3':limit,'pass':max(e1,e2)<.001 and valid and difference<limit})
    report={'scope':'bidirectional nearest-vertex comparison, not a continuous surface Hausdorff proof','input_cache_sha256':source['input_cache_sha256'],'parts':rows,'pass':all(r['pass'] for r in rows)}
    (OUT/'cad_vertex_correspondence.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    if not report['pass']:raise RuntimeError('CAD vertex mismatch')
if __name__=='__main__':run()
