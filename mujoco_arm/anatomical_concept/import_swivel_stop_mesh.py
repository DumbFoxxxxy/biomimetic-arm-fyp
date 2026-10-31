"""Adopt valid native stop CAD as source mesh, then independently recheck it."""
from build_distal_swivels import *
from extract_attachment_regions import projector
from diagnose_swivel_stop_seams import winding
import shutil
def run():
    fid=sys.argv[sys.argv.index('--muscle')+1];folder=DEST/fid;raw=(folder/'canonical_stops.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];v,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);f=index.reshape(-1,3);s=mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)));assert s.status()==mf.Error.NoError and sum(p.volume()>.001 for p in s.decompose())==1
    old,ov,of=load(folder/'distal_socket_base.obj');forward=projector(v,f)(ov)[1];reverse=projector(ov,of)(v)[1];delta=max(float(forward.max()),float(reverse.max()));relative=abs(s.volume()-old.volume())/max(1,old.volume())
    # The rejected MF union retained deep vertices at an internal Boolean seam.
    # Prove they lie inside the valid CAD solid; do not relax an exterior bound.
    seam=np.flatnonzero(forward>=.05);inside=[winding(ov[i],v[f]) for i in seam];assert all(abs(w)>.99 for w in inside),(delta,inside)
    external=max(float(forward[forward<.05].max(initial=0)),float(reverse.max()));symmetric=max(0,float((s-old).volume()))+max(0,float((old-s).volume()));assert external<.05 and relative<1e-4 and symmetric/max(1,old.volume())<1e-4,(external,relative,symmetric)
    archive=folder/'before_native_stop_mesh';archive.mkdir(exist_ok=True)
    for name in ['distal_socket_base.obj','distal_socket_base.stl','cad_input.json']:
        if not (archive/name).exists():shutil.copy2(folder/name,archive/name)
    part=write(folder/'distal_socket_base.stl',s);data=json.loads((folder/'cad_input.json').read_text(encoding='utf8'));data['parts']=[part if p['name']==part['name'] else p for p in data['parts']];(folder/'cad_input.json').write_text(json.dumps(data),encoding='utf8')
    report={'source':'native Boolean union of verified unstopped cap and exact24-sided annular shoulders','source_canonical_STL_sha256':hashlib.sha256(raw).hexdigest(),'relative_volume_change_from_rejected_MF_mesh':relative,'mesh_vertex_surface_deviation_including_invalid_internal_seams_mm':delta,'invalid_internal_seam_vertex_count':len(seam),'internal_seam_vertex_winding_in_valid_CAD_mesh':inside,'exterior_vertex_surface_deviation_mm':external,'material_symmetric_difference_mm3':symmetric,'positive_material_components':1,'OBJ_sha256':hashlib.sha256((folder/'distal_socket_base.obj').read_bytes()).hexdigest(),'fresh_neutral_articulation_assembly_required':True}
    (folder/'canonical_stop_mesh.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(fid,'canonical stop mesh',delta,relative,flush=True)
if __name__=='__main__':run()
