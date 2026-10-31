"""Independent test of numerical bone corridors against actual rigid hardware.

A bone-clear path may hit a terminal cap or its bracket. TPU deformation is
not known, so do not test a moving path against an un-deformed frozen belly.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation
from project_runtime import WorkspacePath as Path

def run():
    assert EXTENSOR_REANCHOR
    source=DEST/'bone_corridor_quick_dense_ports.json';study=json.loads(source.read_text());a=json.loads((DEST/'assembly.json').read_text());j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);rows={r['key']:r for r in a['instances'] if r['kind']!='TPU'};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};inputs=[source,DEST/'assembly.json']+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};records=[]
    for pose in study['poses']:
        q=pose['elbow_deg'];R=Rotation.from_rotvec(axis*np.radians(q)).as_matrix();T=np.c_[R,E-R@E];nodes=np.array(pose['centreline_display_mm']);wire=mf.Manifold.batch_boolean([rod(p,r,.3,24) for p,r in zip(nodes,nodes[1:])]+[mf.Manifold.sphere(.3,16).translate(p) for p in nodes[1:-1]],mf.OpType.Add);hits=[]
        for key,s in world.items():
            target=s.transform(T) if rows[key]['owner']=='elbow' or rows[key].get('nominal_pose_owner')=='elbow' else s
            v=intersection(wire,target)
            if v>.001:hits.append({'part':key,'volume_mm3':v,'nominal_stop_only':rows[key].get('nominal_packing_only',False)})
        links=np.diff(nodes,axis=0);unitlinks=links/np.linalg.norm(links,axis=1)[:,None];turns=np.degrees(np.arccos(np.clip((unitlinks[:-1]*unitlinks[1:]).sum(1),-1,1)))
        records.append({'elbow_deg':q,'rigid_steel_intersections':hits,'max_polyline_turn_deg':float(turns.max()),'proximal_port_turn_deg':float(turns[0]),'distal_port_turn_deg':float(turns[-1]),'clear_all_rigid_envelopes':not hits,'actual_TPU_guide_pose_or_route_verified':False})
        print('corridor hardware',q,'hits',len(hits),'port turn',turns[-1],flush=True)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    (DEST/'bone_corridor_ports_hardware.json').write_text(json.dumps({'scope':'numerical centreline converted to0.6mm steel envelope vs all105 nominal rigid pieces, including two uncrimped stops; actual lattice/guide deformation, friction and drive excluded','source_sha256':hashes,'poses':records,'frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
