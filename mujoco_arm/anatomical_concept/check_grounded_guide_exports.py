"""M7 OBJ/STL roundtrip and actual 109 rigid MuJoCo mesh registration."""
from build_guide_ground_frames import *
from scipy.spatial import cKDTree
import mujoco

def run():
    inputs=json.loads((DEST7/'cad_input.json').read_text());checks=[]
    for r in inputs['parts']:
        name=r['name'];s,v,f=load(DEST7/(name+'.obj'));raw=(DEST7/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,idx=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);restored=mf.Manifold(mf.Mesh(sv.astype(np.float32),idx.reshape(-1,3).astype(np.uint32)));delta=float(abs(tri-v[f]).max());error=abs(restored.volume()-s.volume())/max(1,s.volume());checks.append({'part':name,'coordinate_error_mm':delta,'relative_volume_error':error,'pass':s.status()==restored.status()==mf.Error.NoError and delta<.0001 and error<.0001})
    a=json.loads((DEST7/'assembly.json').read_text());included=set(json.loads((DEST7/'rigid_motion_scope.json').read_text())['included']);m=mujoco.MjModel.from_xml_path(str(DEST7/'arm.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);registration=[]
    for r in a['instances']:
        if r['key'] not in included:continue
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,r['key']);assert gid>=0;mid=m.geom_dataid[gid];start=m.mesh_vertadr[mid];end=start+m.mesh_vertnum[mid];actual=(m.mesh_vert[start:end]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid])*1000;_,v,f=load(ROOT/r['file']);expected=v[np.unique(f)]+np.array(r['translation_display_mm']);error=max(float(cKDTree(actual).query(expected)[0].max()),float(cKDTree(expected).query(actual)[0].max()));registration.append({'part':r['key'],'max_coordinate_error_mm':error,'pass':error<.001})
    assert len(registration)==109 and all(r['pass'] for r in checks+registration)
    result={'scope':'10 new/revised OBJ/STL and109 prescribed rigid neutral model registrations, not installed muscle/drive proof','roundtrips':checks,'registration':registration,'all_pass':True,'assembly_sha256':digest(DEST7/'assembly.json'),'frozen':False};(DEST7/'export_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M7 exports',len(checks),'registered',len(registration),'max mm',max(r['max_coordinate_error_mm'] for r in registration),flush=True)

if __name__=='__main__':run()
