"""Prepare exact audited meshes as design-coordinate mm CAD input."""
import json,hashlib,numpy as np
from verify_joint_candidate import load_candidate,OUT,world
def run():
    report=json.loads((OUT/'retained_shape_rom.json').read_text());cache=OUT/report['mesh_cache']
    if len(report['actions'])!=12 or hashlib.sha256(cache.read_bytes()).hexdigest()!=report['mesh_cache_sha256']:raise RuntimeError('Incomplete/stale audit')
    m,d,s,b=load_candidate(report['mesh_cache']);parts=[]
    for name,solid in s.items():
        actual=world(d,solid,b[name]);mesh=actual.to_mesh();points=np.asarray(mesh.vert_properties)[:,:3];design=points[:,[0,2,1]]*np.array([1,1,-1])
        parts.append({'name':name,'vertices_design_mm':design.tolist(),'faces':np.asarray(mesh.tri_verts).tolist(),'volume_mm3':float(actual.volume())})
    data={'input_cache_sha256':report['mesh_cache_sha256'],'coordinate_system':'mm, X lateral, Y superior, Z posterior; zero assembly placement','parts':parts}
    (OUT/'cad_input.json').write_text(json.dumps(data));print('CAD input prepared',flush=True)
if __name__=='__main__':run()
