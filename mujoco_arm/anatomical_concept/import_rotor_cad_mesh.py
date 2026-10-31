"""Use valid CAD union as the canonical M3 forearm; preserve rejected MF mesh."""
from build_elbow_muscle_mounts import mf,np,load
from build_SC_connection import write
from check_rotor_pocket import DEST
from extract_attachment_regions import projector
import json,hashlib,shutil
def run():
    raw=(DEST/'forearm_from_native.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')])
    triangles=np.frombuffer(raw,dt,offset=84)['v'];v,index=np.unique(triangles.reshape(-1,3),axis=0,return_inverse=True);f=index.reshape(-1,3)
    shape=mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))
    assert shape.status()==mf.Error.NoError,str(shape.status())
    assert len([c for c in shape.decompose() if c.volume()>.001])==1
    old,ov,of=load(DEST/'forearm_with_rotor.obj')
    forward=projector(v,f)(ov)[1];reverse=projector(ov,of)(v)[1]
    delta=max(float(forward.max()),float(reverse.max()));relative=abs(shape.volume()-old.volume())/old.volume()
    assert delta<.1 and relative<1e-4,(delta,relative)
    archive=DEST/'before_native_canonical_mesh';archive.mkdir(exist_ok=True)
    for name in ['forearm_with_rotor.obj','forearm_with_rotor.stl','cad_input.json','bone_motion_audit.json','selection.json','native_union_diagnostic.json']:
        if (DEST/name).exists() and not (archive/name).exists():shutil.copy2(DEST/name,archive/name)
    part=write(DEST/'forearm_with_rotor.stl',shape)
    data=json.loads((DEST/'cad_input.json').read_text(encoding='utf8'));data['parts']=[part if p['name']=='forearm_with_rotor' else p for p in data['parts']]
    (DEST/'cad_input.json').write_text(json.dumps(data),encoding='utf8')
    report={'reason':'MF Boolean union mesh failed sewn BRep; native24-sided prismatic rotor fused with verified parent STEP is valid. Canonical mesh now re-exported from this CAD union, no tolerance relaxation.','source_native_STL_sha256':hashlib.sha256(raw).hexdigest(),'source_STEP_sha256':hashlib.sha256((DEST/'cad/forearm_with_rotor.step').read_bytes()).hexdigest(),'mesh_vertex_to_other_mesh_surface_max_mm':delta,'relative_volume_change_from_rejected_mesh':relative,'positive_material_components':1,'fresh_motion_and_neutral_checks_required':True,'frozen':False}
    (DEST/'canonical_mesh_revision.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)
if __name__=='__main__':run()
