"""Independent export topology and neutral soft/reference-bone checks."""
from build_direct_muscle_yoke import *
from scipy.spatial import cKDTree
import mujoco

def run():
    assert '--access-relief' in sys.argv or EXTENSOR_REANCHOR
    a=json.loads((DEST/'assembly.json').read_text());inputs=json.loads((DEST/'cad_input.json').read_text());parts=[]
    for p in inputs['parts']:
        # Original five mesh pairs remain at their recorded original paths.
        folder=DEST if EXTENSOR_REANCHOR or p['name']=='direct_muscle_yoke' else OUT/'direct_muscle_interfaces/mechanical_yoke_shifted'
        s,v,f=load(folder/(p['name']+'.obj'));raw=(folder/(p['name']+'.stl')).read_bytes()
        dtype=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dtype,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True)
        stl=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)))
        error=abs(stl.volume()-s.volume())/max(1,abs(s.volume()));vertex_error=max(cKDTree(v).query(sv)[0].max(),cKDTree(sv).query(v)[0].max())
        passed=s.status()==stl.status()==mf.Error.NoError and error<1e-4 and vertex_error<.001
        parts.append({'part':p['name'],'OBJ_sha256':hashlib.sha256((folder/(p['name']+'.obj')).read_bytes()).hexdigest(),'STL_sha256':hashlib.sha256(raw).hexdigest(),'STL_relative_volume_error':float(error),'bidirectional_vertex_error_mm':float(vertex_error),'pass':bool(passed)})
    refs=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'));bones=[]
    for r in refs['models']:
        if r['kind']=='bone_reference':bones.append((r['id'],load(OUT/'repaired_supports'/(r['id']+'.obj'),raw=True)[0]))
    hits=[];checked=[]
    for r in a['instances']:
        if r['kind']!='TPU' and not r.get('nominal_packing_only'):continue
        s=load(ROOT/r['file'])[0].translate(r['translation_display_mm']);checked.append(r['key'])
        for name,bone in bones:
            volume=intersection(s,bone)
            if volume>.001:hits.append({'pair':[r['key'],name],'volume_mm3':volume,'classification':'neutral_soft_review' if r['kind']=='TPU' else 'nominal_stop_review'})
    registrations={}
    for name in ['neutral_view.xml','arm.xml']:
        m=mujoco.MjModel.from_xml_path(str(DEST/name));d=mujoco.MjData(m);mujoco.mj_forward(m,d);scene=[]
        for r in a['instances']:
            gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,r['key'])
            if gid<0:
                assert name=='arm.xml' and (r['kind']=='TPU' or r.get('nominal_packing_only'))
                continue
            mid=m.geom_dataid[gid];start=m.mesh_vertadr[mid];end=start+m.mesh_vertnum[mid]
            actual=(m.mesh_vert[start:end]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid])*1000
            _,v,f=load(ROOT/r['file']);expected=v[np.unique(f)]+np.array(r['translation_display_mm'])
            err=max(cKDTree(actual).query(expected)[0].max(),cKDTree(expected).query(actual)[0].max())
            scene.append({'part':r['key'],'error_mm':float(err),'pass':bool(err<.001)})
        registrations[name]=scene
    report={'scope':'six generated rigid OBJ/STL pairs, two TPU and two nominal stops against reference bones at neutral, and compiled MuJoCo geometry registration; no deformation or motion claim','parts':parts,'all_roundtrips_pass':all(p['pass'] for p in parts),'reference_bone_count':len(bones),'neutral_checked_instances':checked,'neutral_reference_intersections':hits,'compiled_scenes':registrations,'frozen':False}
    (DEST/'mesh_export_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('six rigid mesh roundtrips',report['all_roundtrips_pass'],'neutral refs',len(bones),'hits',len(hits),flush=True)

if __name__=='__main__':run()
