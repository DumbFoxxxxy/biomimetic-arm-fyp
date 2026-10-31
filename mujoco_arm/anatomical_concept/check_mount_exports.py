"""M2 STL/OBJ and compiled neutral/rigid scene registration checks."""
from build_elbow_muscle_mounts import *
from scipy.spatial import cKDTree
import mujoco

def run():
    assembly=json.loads((DEST/'assembly.json').read_text());checks=[]
    for name in assembly['generated_components']:
        v,f=read_obj(str(DEST/(name+'.obj')));source=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)));raw=(DEST/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);re=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));relative=abs(re.volume()-source.volume())/max(1,abs(source.volume()));error=float(np.max(abs(tri-v[f])))
        checks.append({'part':name,'OBJ_status':str(source.status()),'STL_status':str(re.status()),'relative_volume_error':relative,'coordinate_error_mm':error,'pass':source.status()==re.status()==mf.Error.NoError and relative<1e-4 and error<1e-4})
    scenes={}
    scope=json.loads((DEST/'rigid_motion_scope.json').read_text());rigid=set(scope['included'])
    for filename,selected in [('neutral_view.xml',assembly['instances']),('arm.xml',[r for r in assembly['instances'] if r['key'] in rigid])]:
        m=mujoco.MjModel.from_xml_path(str(DEST/filename));d=mujoco.MjData(m);mujoco.mj_forward(m,d);rows=[]
        for r in selected:
            gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,r['key'])
            if gid<0:raise RuntimeError('Missing scene geometry '+r['key'])
            mid=int(m.geom_dataid[gid]);start=int(m.mesh_vertadr[mid]);n=int(m.mesh_vertnum[mid]);compiled=m.mesh_vert[start:start+n]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid]
            v,f=read_obj(str(ROOT/r['file']));expected=(v[np.unique(f)]+np.array(r['translation_display_mm']))/1000
            delta=max(float(cKDTree(expected).query(compiled)[0].max()),float(cKDTree(compiled).query(expected)[0].max()))*1000
            rows.append({'part':r['key'],'max_mesh_coordinate_error_mm':delta,'pass':delta<.001})
        scenes[filename]=rows
    result={'scope':'export fidelity and zero-pose registration, not force/compliance/motion proof','parts':checks,'scenes':scenes,'all_generated_roundtrips_pass':all(r['pass'] for r in checks),'all_scene_registration_pass':all(r['pass'] for rows in scenes.values() for r in rows),'frozen':False}
    (DEST/'export_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    if not result['all_generated_roundtrips_pass'] or not result['all_scene_registration_pass']:raise RuntimeError('M2 export mismatch')
    print('M2 export parts',len(checks),'scene counts',{k:len(v) for k,v in scenes.items()},'max mm',max(r['max_mesh_coordinate_error_mm'] for rows in scenes.values() for r in rows))

if __name__=='__main__':run()
