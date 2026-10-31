"""Exact M5 OBJ/STL fidelity, positive connectivity and131-scene registration."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_distal_swivels import *
from scipy.spatial import cKDTree
import mujoco
from scipy.spatial.transform import Rotation
def run():
    checks=[]
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=json.loads((folder/'design.json').read_text(encoding='utf8'))
        for name in design['parts']:
            s,v,f=load(folder/(name+'.obj'));raw=(folder/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);re=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));relative=abs(re.volume()-s.volume())/max(1,abs(s.volume()));error=float(abs(tri-v[f]).max());positive=sum(c.volume()>.001 for c in re.decompose())
            check={'muscle':fid,'part':name,'STL_status':str(re.status()),'relative_volume_error':relative,'coordinate_error_mm':error,'positive_material_components':positive,'OBJ_sha256':hashlib.sha256((folder/(name+'.obj')).read_bytes()).hexdigest(),'STL_sha256':hashlib.sha256(raw).hexdigest(),'pass':re.status()==mf.Error.NoError and relative<1e-4 and error<1e-4 and positive==1};checks.append(check)
        C=np.array(design['selected']['ball_centre_display_mm']);u=np.array(design['selected']['neutral_terminal_axis']);wire=rod(C+u*3.5,C+u*9.5,.3);terminal=load(folder/'ball_crimp_blank.obj')[0]
        assert intersection(wire,terminal)<=.001
        minimum=float(wire.min_gap(terminal,3));(folder/'wire_stub_check.json').write_text(json.dumps({'scope':'straight0.6mm cable inside uncrimped terminal bore, no external route or grip claim','intersection_mm3':intersection(wire,terminal),'minimum_gap_mm_capped3':minimum,'actual_crimp_grip_verified':False},indent=2),encoding='utf8')
    assembly=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));m=_load_project_model(str(DEST/'neutral_view.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);scene=[]
    for r in assembly['instances']:
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,r['key']);assert gid>=0;mid=int(m.geom_dataid[gid]);start=int(m.mesh_vertadr[mid]);count=int(m.mesh_vertnum[mid]);actual=m.mesh_vert[start:start+count]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid];_,v,f=load(ROOT/r['file']);expected=(v[np.unique(f)]+r['translation_display_mm'])/1000
        delta=max(float(cKDTree(expected).query(actual)[0].max()),float(cKDTree(actual).query(expected)[0].max()))*1000;scene.append({'part':r['key'],'mesh_registration_error_mm':delta,'pass':delta<.001})
    animation=[]
    for fid in ['FJ1486','FJ1477']:
        design=json.loads((DEST/fid/'design.json').read_text(encoding='utf8'));p=design['selected'];C=np.array(p['ball_centre_display_mm']);u=np.array(p['neutral_terminal_axis']);v=np.array(p['retainer_v_axis']);source,vertices,faces=load(DEST/fid/'ball_crimp_blank.obj')
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,fid+'_swivel_ball_crimp_blank');bid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,fid+'_swivel_ball_crimp_blank_motion');mocap=int(m.body_mocapid[bid]);mid=int(m.geom_dataid[gid]);start=int(m.mesh_vertadr[mid]);count=int(m.mesh_vertnum[mid])
        for label,axis in [('v',v),('w',np.cross(u,v))]:
            for angle in [-25,0,25]:
                d.mocap_pos[:]=0;d.mocap_quat[:]=[1,0,0,0];rotation=Rotation.from_rotvec(axis*np.radians(angle));R=rotation.as_matrix();quat=rotation.as_quat();d.mocap_pos[mocap]=(C-R@C)/1000;d.mocap_quat[mocap]=quat[[3,0,1,2]];mujoco.mj_forward(m,d)
                actual=m.mesh_vert[start:start+count]@d.geom_xmat[gid].reshape(3,3).T+d.geom_xpos[gid];expected=((vertices-C)@R.T+C)/1000
                delta=max(float(cKDTree(expected).query(actual)[0].max()),float(cKDTree(actual).query(expected)[0].max()))*1000;animation.append({'muscle':fid,'axis':label,'angle_deg':angle,'rotation_registration_error_mm':delta,'pass':delta<.001})
    assert len(checks)==22 and len(scene)==131 and all(r['pass'] for r in checks+scene+animation)
    report={'scope':'manufacturing mesh fidelity, neutral scene registration and prescribed terminal rotation registration; no drive/deformation/force proof','parts':checks,'scene':scene,'prescribed_terminal_pose_checks':animation,'neutral_scene_sha256':hashlib.sha256((DEST/'neutral_view.xml').read_bytes()).hexdigest(),'pass':True,'frozen':False};(DEST/'export_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('M5 exports22 scene131, max mm',max(r['mesh_registration_error_mm'] for r in scene),'prescribed terminal poses',len(animation))
if __name__=='__main__':run()
