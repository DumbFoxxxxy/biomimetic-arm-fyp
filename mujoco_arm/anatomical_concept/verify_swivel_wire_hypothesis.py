"""Independent endpoint, clearance and local self-contact checks of wire studies.

Verifies hypotheses against actual files; never promotes them to full TPU or
actuation validation. Failed cases are retained as evidence, not assertions of
an impossible muscle design.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from project_runtime import WorkspacePath as Path

def cable(points,radius=.3):
    return mf.Manifold.batch_boolean([rod(p,q,radius,24) for p,q in zip(points,points[1:])]+[mf.Manifold.sphere(radius,16).translate(p) for p in points[1:-1]],mf.OpType.Add)

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';candidate=parent/'captive_keepers';folder=candidate/('wire_pose_regularized_fine' if '--fine' in sys.argv else 'wire_pose_regularized_critical' if '--regularized' in sys.argv else 'wire_pose_study');path=folder/'wire_guides_quick.json';study=json.loads(path.read_text());study_digest=hashlib.sha256(path.read_bytes()).hexdigest();a=json.loads((candidate/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};motions=json.loads((parent/'conditional_motion.json').read_text())['poses'];world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if r['kind']!='TPU'}
    for p,h in study['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,(p,'stale source')
    tubes={};data={};records=[]
    for fid in set(p['muscle'] for p in study['poses']):
        d=json.loads((muscle_folder(fid)/'design.json').read_text());offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));body=load(ROOT/rows[fid+'_TPU_muscle']['file'])[0].translate(rows[fid+'_TPU_muscle']['translation_display_mm']);original=[]
        for g in d['actual_short_tube_guides']:
            c=np.array(g['centre_display_mm'])+offset;t=np.array(g['axis_display']);half=g['length_mm']/2;original.append(body^rod(c-t*half,c+t*half,g['outer_radius_mm']+.015,48))
        tubes[fid]=original;data[fid]=(d,offset)
    for case in study['poses']:
        fid=case['muscle'];q=case['elbow_deg'];pose=next(p for p in motions if p['elbow_deg']==q);Tc=np.array(pose['cap_poses'][fid]);Tf=np.array(pose['forearm_pose']);d,offset=data[fid];A,B=np.array(d['pad_centres_display_mm'])+offset;z=np.array(d['basis_columns_display'])[:,2];P=np.array(case['centreline_display_mm']);expected=np.array([A-z*5.15,A-z*2.15,*((np.array([B+z*2.15,B+z*4.6])@Tc[:,:3].T+Tc[:,3]).tolist())]);error=float(abs(P[[0,1,-2,-1]]-expected).max());assert error<1e-7
        wire=cable(P);guides=[];tube_wire=[];gap_guide=2.
        for i,item in enumerate(case['guide_checks']):
            T=np.array(item['pose_affine_mm']);assert np.linalg.norm(T[:,:3].T@T[:,:3]-np.eye(3))<1e-8 and abs(np.linalg.det(T[:,:3])-1)<1e-8;centre=np.array(d['actual_short_tube_guides'][i]['centre_display_mm'])+offset;assert np.linalg.norm(T[:,:3]@centre+T[:,3]-P[i+2])<1e-7;g=tubes[fid][i].transform(T);guides.append(g);volume=intersection(wire,g)
            if volume>.001:tube_wire.append({'guide':i,'volume_mm3':volume})
            gap_guide=min(gap_guide,float(wire.min_gap(g,2)))
        guide_self=[]
        for i,j in itertools.combinations(range(len(guides)),2):
            v=intersection(guides[i],guides[j])
            if v>.001:guide_self.append({'guides':[i,j],'volume_mm3':v})
        segments=[cable([p,r]) for p,r in zip(P,P[1:])];wire_self=[]
        for i,j in itertools.combinations(range(len(segments)),2):
            if j<=i+1:continue
            v=intersection(segments[i],segments[j])
            if v>.001:wire_self.append({'segments':[i,j],'volume_mm3':v})
        hard=[];gap_hard=2.
        for key,body in world.items():
            owner=rows[key]['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;s=body.transform(T) if T is not None else body;v=intersection(wire,s)
            if v>.001:hard.append({'part':key,'volume_mm3':v})
            aa=np.array(wire.bounding_box());bb=np.array(s.bounding_box());lower=np.linalg.norm(np.maximum(0,np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])))
            if lower<2:gap_hard=min(gap_hard,float(wire.min_gap(s,2)))
        record={'muscle':fid,'elbow_deg':q,'port_endpoint_error_mm':error,'steel_rigid_intersections':hard,'steel_guide_material_intersections':tube_wire,'guide_piece_self_intersections':guide_self,'wire_nonadjacent_self_intersections':wire_self,'minimum_steel_to_guide_gap_mm_capped2':gap_guide,'minimum_steel_to_rigid_gap_mm_capped2':gap_hard,'complete_TPU_mesh_deformation_verified':False,'continuous_driverange_verified':False};records.append(record);print(fid,q,'hard',len(hard),'guide',len(tube_wire),'self',len(guide_self),len(wire_self),'guidegap',gap_guide,flush=True)
    assert hashlib.sha256(path.read_bytes()).hexdigest()==study_digest
    for p,h in study['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    (folder/'independent_review.json').write_text(json.dumps({'scope':'actual endpoint+rigid metal+manufactured short-guide lumen andlocal self-contact checks; no complete TPU equilibrium, body contact,friction or actuation proof; cases do not establish a continuous range','source_study_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_sha256':study['source_sha256'],'poses':records,'frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
