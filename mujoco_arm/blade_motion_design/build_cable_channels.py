"""Subtract finite-radius swept cable corridors from the repaired scapula.

Channels are free-clearance corridors; no force-bearing wall guidance is
assumed. Joint-interface regions are protected and checked explicitly.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys,json,argparse
import numpy as np
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parent.parent
SOURCE=PROJECT/'tmp/shoulder_motion_design_20261029';EXPORT=ROOT.parent/'exports/blade_motion_design'
sys.path[:0]=[str(SOURCE/'vendor'),str(ROOT.parent/'cable_design_v8'),str(ROOT.parent)]
import manifold3d as md
import v8_shoulder as v8
from stl_io import write_binary_stl
sys.stdout.reconfigure(encoding='utf-8')

def cylinder_between(a,b,radius=2.):
    direction=np.asarray(b)-a;length=float(np.linalg.norm(direction));z=direction/length
    helper=np.array([1.,0,0]) if abs(z[0])<.8 else np.array([0.,1,0])
    x=np.cross(helper,z);x/=np.linalg.norm(x);y=np.cross(z,x)
    return md.Manifold.cylinder(length,radius,circular_segments=24).transform(np.column_stack([x,y,z,a]))

def clip_segment(a,b,lo,hi):
    d=b-a;t0=0.;t1=1.
    for axis in range(3):
        if abs(d[axis])<1e-12:
            if a[axis]<lo[axis] or a[axis]>hi[axis]:return None
        else:
            u=(lo[axis]-a[axis])/d[axis];v=(hi[axis]-a[axis])/d[axis]
            t0=max(t0,min(u,v));t1=min(t1,max(u,v))
            if t0>=t1:return None
    return a+t0*d,a+t1*d

def save_solid(path,solid):
    mesh=solid.to_mesh();points=np.asarray(mesh.vert_properties)[:,:3];faces=np.asarray(mesh.tri_verts)
    mapped=points@v8.DESIGN_TO_MJ.T*.001
    write_binary_stl(path,mapped[faces],'Blade shoulder: real swept cable corridor subtraction; MJ local metres')

def build_channels(hybrid=False):
    cfg=v8.load_config();fit=json.loads((SOURCE/'paper_motion_fit_repaired.json').read_text());anchors=json.loads((SOURCE/('anchor_search_hybrid_channels.json' if hybrid else 'anchor_search_no_contact_repaired_channels.json')).read_text())
    if not anchors['pass']:raise RuntimeError('Channel candidate does not pass pull-only directional test')
    model=v8.compile_tree(v8.build_shoulder_tree(cfg,0));gid=v8.object_id(model,v8.mujoco.mjtObj.mjOBJ_GEOM,'scapula_mesh');mid=int(model.geom_dataid[gid])
    verts=v8.scapula_mesh_vertices_design(cfg)+12*v8.source_basis(cfg)[2]
    fa=int(model.mesh_faceadr[mid]);fn=int(model.mesh_facenum[mid]);faces=model.mesh_face[fa:fa+fn].copy()
    solid=md.Manifold(md.Mesh(np.ascontiguousarray(verts,dtype=np.float32),np.ascontiguousarray(faces,dtype=np.uint32)))
    print('SOURCE',solid.status(),solid.volume(),len(solid.decompose()),flush=True)
    if solid.is_empty():raise RuntimeError('Source scapula is not a valid closed manifold')
    zero=np.array(fit['poses'][0]['q_rad']);q_samples=[]
    for name in ['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']:
        seq=[zero]+[np.array(p['q_rad']) for p in fit['poses'] if p.get('action')==name]
        for a,b in zip(seq[:-1],seq[1:]):
            for f in np.linspace(0.,1.,5,endpoint=False):q_samples.append(a*(1-f)+b*f)
        q_samples.append(seq[-1])
    fks=[v8.direct_kinematics(cfg,q) for q in q_samples]
    radius=2.;lo=verts.min(0)-radius-1.;hi=verts.max(0)+radius+1.
    corridors=[];route_info=[]
    for row in anchors['routes']:
        if hybrid and row['id'] not in anchors['channel_cables']:continue
        cuts=[];end=np.array(row['eye_local_ac_mm']);fixed=np.array(row['fixed_eye_world_mm'])
        for fk in fks:
            start=fk['body_rotation'].T@(fixed-fk['ac_mm']);clipped=clip_segment(start,end,lo,hi)
            if clipped is not None:cuts.append(cylinder_between(*clipped,radius))
        corridor=md.Manifold.batch_boolean(cuts,md.OpType.Add)
        actual=md.Manifold.batch_boolean([solid,corridor],md.OpType.Intersect)
        corridors.append(corridor)
        info={'cable':row['id'],'sweep_samples':len(fks),'cutter_segments':len(cuts),'removed_candidate_volume_mm3':actual.volume(),
              'corridor_radius_mm':radius,'actual_cut_bounds_ac_local_mm':list(actual.bounding_box()) if not actual.is_empty() else None}
        route_info.append(info);print('CHANNEL',json.dumps(info),flush=True)
    union=md.Manifold.batch_boolean(corridors,md.OpType.Add)
    ac=np.array(cfg['geometry_mm']['AC']);gh=np.array([159.06,187.54,7.30])-ac;w=v8.source_basis(cfg)[2]
    protect=md.Manifold.batch_boolean([md.Manifold.sphere(10.,circular_segments=32).translate(center) for center in [np.zeros(3),gh,12*w,gh+12*w]],md.OpType.Add)
    relevant=md.Manifold.batch_boolean([solid,union],md.OpType.Intersect)
    conflict=md.Manifold.batch_boolean([relevant,protect],md.OpType.Intersect)
    allowed=md.Manifold.batch_boolean([union,protect],md.OpType.Subtract)
    result=md.Manifold.batch_boolean([solid,allowed],md.OpType.Subtract)
    EXPORT.mkdir(parents=True,exist_ok=True)
    prefix='scapula_hybrid' if hybrid else 'scapula_channels'
    save_solid(EXPORT/f'{prefix}_body_local.stl',result)
    save_solid(EXPORT/'scapula_uncut_repaired_body_local.stl',solid)
    save_solid(EXPORT/f'{prefix}_removal_body_local.stl',md.Manifold.batch_boolean([solid,result],md.OpType.Subtract))
    mesh=result.to_mesh();points=np.asarray(mesh.vert_properties)[:,:3];tri=np.asarray(mesh.tri_verts)
    np.savez_compressed(EXPORT/f'{prefix}_design_geometry.npz',vertices_ac_local_mm=points,faces=tri)
    report={'scope':'Actual manifold boolean subtraction of radius-2mm wire sweeps; preserves radius-10mm joint zones',
        'source_status':str(solid.status()),'result_status':str(result.status()),'source_volume_mm3':solid.volume(),'result_volume_mm3':result.volume(),
        'removed_volume_mm3':solid.volume()-result.volume(),'removed_fraction':1-result.volume()/solid.volume(),
        'source_components':len(solid.decompose()),'result_components':len(result.decompose()),
        'protected_joint_conflict_mm3':conflict.volume(),'channel_routes':route_info,
        'geometry_boolean_pass':bool(not result.is_empty() and len(result.decompose())==len(solid.decompose())),
        'free_cable_clearance_verified':False,'ready_for_manufacture':False}
    (EXPORT/('hybrid_channel_boolean_report.json' if hybrid else 'channel_boolean_report.json')).write_text(json.dumps(report,indent=2));print('SUMMARY',json.dumps({k:v for k,v in report.items() if k!='channel_routes'}),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--hybrid',action='store_true');args=parser.parse_args();build_channels(args.hybrid)
