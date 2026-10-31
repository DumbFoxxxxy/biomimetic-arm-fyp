"""Actual mesh and MuJoCo pose registration, independent of collision sampling."""
from build_direct_muscle_yoke import *
from scipy.spatial import cKDTree
import mujoco
from pathlib import Path

def run():
    dest=DEST/'paired_swivel_support'
    if '--keepers' in sys.argv:dest=dest/'captive_keepers'
    a=json.loads((dest/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};cad=json.loads((dest/'cad_input.json').read_text());exports=[]
    for p in cad['parts']:
        path=ROOT/rows[p['name']]['file'];s,v,f=load(path);raw=path.with_suffix('.stl').read_bytes();dtype=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dtype,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);stl=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));err=abs(stl.volume()-s.volume())/max(1,abs(s.volume()));distance=max(cKDTree(v).query(sv)[0].max(),cKDTree(sv).query(v)[0].max());passed=s.status()==stl.status()==mf.Error.NoError and err<1e-4 and distance<.001;assert passed,p['name'];exports.append({'part':p['name'],'STL_relative_volume_error':float(err),'bidirectional_vertex_error_mm':float(distance),'pass':bool(passed)})
    scenes={};motion=json.loads(((dest.parent if '--keepers' in sys.argv else dest)/'conditional_motion.json').read_text())['poses'];source={k:load(ROOT/r['file'])[1]+np.array(r['translation_display_mm']) for k,r in rows.items()}
    for scene in ['neutral_view.xml','rigid_hinges.xml']:
        m=mujoco.MjModel.from_xml_path(str(dest/scene));d=mujoco.MjData(m);cases=[('nominal_zero',None)]
        if scene=='rigid_hinges.xml':cases += [(f'conditional_{q}',next(p for p in motion if p['elbow_deg']==q)) for q in [0,90,150]]
        scene_results=[]
        for label,pose in cases:
            d.qpos[:]=0
            if pose:
                for name,value in [('elbow',pose['elbow_deg'])]+[(fid+'_swivel',pose['relative_cap_angles_deg'][fid]) for fid in ['FJ1486','FJ1477']]:
                    jid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name);assert jid>=0;d.qpos[m.jnt_qposadr[jid]]=np.radians(value)
            mujoco.mj_forward(m,d);results=[]
            for key,r in rows.items():
                gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,key)
                if gid<0:assert scene=='rigid_hinges.xml' and (r['kind']=='TPU' or r.get('nominal_packing_only'));continue
                mid=m.geom_dataid[gid];start=m.mesh_vertadr[mid];end=start+m.mesh_vertnum[mid];actual=(m.mesh_vert[start:end]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid])*1000;expected=source[key]
                if pose:
                    owner=r['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else np.array(pose['forearm_pose']) if owner=='elbow' else None
                    if T is not None:expected=expected@T[:,:3].T+T[:,3]
                error=max(cKDTree(actual).query(expected)[0].max(),cKDTree(expected).query(actual)[0].max());assert error<.001,(scene,label,key,error);results.append({'part':key,'error_mm':float(error),'pass':True})
            scene_results.append({'pose':label,'registered_count':len(results),'maximum_error_mm':max(r['error_mm'] for r in results),'parts':results})
        scenes[scene]=scene_results
    (dest/'export_checks.json').write_text(json.dumps({'scope':f'{len(exports)}manufactured rigid mesh roundtrips andactual MuJoCo source-coordinate registration atnominal andthreeconditional poses; no drive/contact physics validation','parts':exports,'compiled_scenes':scenes,'assembly_sha256':hashlib.sha256((dest/'assembly.json').read_bytes()).hexdigest(),'frozen':False},indent=2),encoding='utf8');print('paired mesh export and all compiled scene poses registered',flush=True)

if __name__=='__main__':run()
