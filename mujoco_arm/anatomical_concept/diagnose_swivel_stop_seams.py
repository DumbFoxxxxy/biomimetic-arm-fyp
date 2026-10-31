"""Distinguish invalid internal Boolean seams from actual exterior changes."""
from build_distal_swivels import *
from extract_attachment_regions import projector
def winding(p,tri):
    a=tri[:,0]-p;b=tri[:,1]-p;c=tri[:,2]-p;la=np.linalg.norm(a,axis=1);lb=np.linalg.norm(b,axis=1);lc=np.linalg.norm(c,axis=1);num=np.einsum('ij,ij->i',a,np.cross(b,c));den=la*lb*lc+np.einsum('ij,ij->i',a,b)*lc+np.einsum('ij,ij->i',b,c)*la+np.einsum('ij,ij->i',c,a)*lb
    return float((2*np.arctan2(num,den)).sum()/(4*np.pi))
def run():
    folder=DEST/'FJ1486';raw=(folder/'canonical_stops.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];v,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);f=index.reshape(-1,3);s=mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)));old,ov,of=load(folder/'distal_socket_base.obj');q,d,indices=projector(v,f)(ov);reverse=projector(ov,of)(v)[1];bad=np.flatnonzero(d>=.05)
    report={'old_decomposed_volumes_mm3':[float(c.volume()) for c in old.decompose()],'forward_max_mm':float(d.max()),'reverse_max_mm':float(reverse.max()),'new_minus_old_mm3':float((s-old).volume()),'old_minus_new_mm3':float((old-s).volume()),'bad_vertices':[{'vertex':int(i),'point_mm':ov[i].tolist(),'nearest_new_surface_mm':q[i].tolist(),'distance_mm':float(d[i]),'winding_in_new_solid':winding(ov[i],v[f])} for i in bad]}
    (folder/'stop_seam_diagnostic.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)
if __name__=='__main__':run()
