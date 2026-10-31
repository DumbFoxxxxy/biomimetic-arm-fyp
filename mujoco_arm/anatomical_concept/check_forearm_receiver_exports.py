"""M8 manufacturing roundtrip and actual compiled neutral receiver meshes."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_forearm_receivers import *
from scipy.spatial import cKDTree
import mujoco
MOUNT=DEST/'rephased_anchors'

def run():
    selected=[]
    for fid in ['FJ1486','FJ1477']:
        for r in json.loads((DEST/fid/'cad_input.json').read_text())['parts']:selected.append(DEST/fid/r['name'])
    selected.append(MOUNT/'lever_receiving_spotfaces');roundtrips=[]
    for path in selected:
        solid,v,f=load(path.with_suffix('.obj'));raw=path.with_suffix('.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,idx=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);s=mf.Manifold(mf.Mesh(sv.astype(np.float32),idx.reshape(-1,3).astype(np.uint32)));error=float(abs(tri-v[f]).max());relative=abs(s.volume()-solid.volume())/max(1,solid.volume());roundtrips.append({'part':str(path.relative_to(DEST)),'coordinate_error_mm':error,'relative_volume_error':relative,'pass':s.status()==solid.status()==mf.Error.NoError and error<.0001 and relative<.0001})
    xml=MOUNT/'whole_arm_direction_study.xml';m=_load_project_model(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d);scene={'side_lever':MOUNT/'lever_receiving_spotfaces.obj'}
    for path in MOUNT.glob('*_neutral.obj'):scene[path.stem]=path
    registrations=[]
    for key,path in scene.items():
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,key);assert g>=0;mi=m.geom_dataid[g];start=m.mesh_vertadr[mi];end=start+m.mesh_vertnum[mi];actual=(m.mesh_vert[start:end]@d.geom_xmat[g].reshape(3,3).T+d.geom_xpos[g])*1000;_,v,f=load(path);expected=v[np.unique(f)];error=max(float(cKDTree(actual).query(expected)[0].max()),float(cKDTree(expected).query(actual)[0].max()));registrations.append({'part':key,'max_coordinate_error_mm':error,'pass':error<.001})
    assert len(roundtrips)==len(registrations)==27 and all(r['pass'] for r in roundtrips+registrations)
    (MOUNT/'export_checks.json').write_text(json.dumps({'scope':'27 new/revised receiver manufacturing meshes and actual static MuJoCo registration only, no drive proof','roundtrips':roundtrips,'registration':registrations,'all_pass':True,'source_xml_sha256':hashlib.sha256(xml.read_bytes()).hexdigest(),'frozen':False},indent=2),encoding='utf8');print('M8 roundtrips27,registered27,max mm',max(r['max_coordinate_error_mm'] for r in registrations),flush=True)

if __name__=='__main__':run()
