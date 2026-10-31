"""Print-only placements, preserving the separate common-coordinate assembly files."""
from parameterize_elbow_muscles import *
from build_SC_connection import write

def run(muscles=('FJ1486','FJ1477')):
    for fid in muscles:
        folder=DEST/fid;report=json.loads((folder/'design.json').read_text(encoding='utf8'));F=np.array(report['basis_columns_display']);E=np.array(report['coordinate_origin_display_mm']);out=folder/'print_candidates';out.mkdir(exist_ok=True);rows=[]
        for name in ['TPU_muscle','proximal_rigid_cap','distal_rigid_cap']:
            v,f=read_obj(str(folder/(name+'.obj')));s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            P=np.array([[1,0,0],[0,0,-1],[0,1,0]]) if name=='TPU_muscle' else np.diag([1,-1,-1]) if name=='proximal_rigid_cap' else np.eye(3)
            R=P@F.T;local=(v-E)@R.T;shift=np.r_[-(local.max(0)[:2]+local.min(0)[:2])/2,-local[:,2].min()];translation=shift-R@E
            target=out/(name+'.stl');placed=s.transform(np.c_[R,translation]);write(target,placed)
            raw=target.read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];pv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);reload=mf.Manifold(mf.Mesh(pv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));relative=abs(reload.volume()-s.volume())/max(1,abs(s.volume()))
            passed=reload.status()==mf.Error.NoError and relative<1e-4
            if not passed:raise RuntimeError('Print STL roundtrip failed '+fid+' '+name)
            rows.append({'part':name,'source_OBJ_sha256':hashlib.sha256((folder/(name+'.obj')).read_bytes()).hexdigest(),'print_STL_sha256':hashlib.sha256(raw).hexdigest(),'print_STL_roundtrip_pass':passed,'print_STL_relative_volume_error':relative,'original_to_print_affine_mm':np.c_[R,translation].tolist(),'volume_change_mm3':float(placed.volume()-s.volume()),'print_bounds_mm':list(placed.bounding_box()),'orientation':'20mm nonload-bearing cell axis vertical, contraction axis horizontal; candidate only, support/brim needs slicer verification' if name=='TPU_muscle' else 'flat TPU contact face toward bed; test bore and eyelet quality','material':'TPU, paper reference95A; actual batch not characterized' if name=='TPU_muscle' else 'rigid material not yet selected/validated'})
        data={'scope':'print placement candidates only; no slicing, support, tolerance, strength or printer qualification. Do not use these shifted files for anatomical assembly.','units':'mm','parts':rows,'fasteners':'metal reference envelopes remain in parent directory; print set contains only TPU body and rigid caps','frozen':False}
        (out/'print_transforms.json').write_text(json.dumps(data,indent=2),encoding='utf8');print(fid,'print placements',len(rows))

if __name__=='__main__':run()
