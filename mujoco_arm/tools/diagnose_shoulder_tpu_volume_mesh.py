"""Mesher control test and localized diagnostics of complete-muscle failures.

A successful cube control only validates the import/volume API path. It
does not promote failed complete TPU models to elasticity-ready status.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
import build_shoulder_tpu_volume_mesh as M

BASE=A.DESIGN/'s1_y2_transition_stock/volume_mesh'

def segment_facet(p,q,triangle):
    a,b,c=triangle;direction=q-p;e1=b-a;e2=c-a
    cross=np.cross(direction,e2);det=float(e1@cross)
    scale=np.linalg.norm(direction)*np.linalg.norm(e1)*np.linalg.norm(e2)
    if abs(det)<=max(1e-30,scale*1e-12):return None
    offset=p-a;u=float(offset@cross/det);v=float(direction@np.cross(offset,e1)/det)
    t=float(e2@np.cross(offset,e1)/det)
    # Strict facet interior distinguishes a crossing/T-junction from intended
    # shared triangle edges; segment endpoints are permitted as PLC defects.
    if 1e-7<u and 1e-7<v and u+v<1-1e-7 and -1e-7<=t<=1+1e-7:
        return dict(t=t,u=u,v=v,point_local_mm=(p+direction*t).tolist())
    return None

def localize(folder,point):
    data=np.load(folder/'input_boundary.npz');V=data['nodes_local_mm'];F=data['triangles'];tri=V[F]
    lower=tri.min(axis=1);upper=tri.max(axis=1);radius=.2
    selected=np.flatnonzero(np.all(lower<=point+radius,axis=1)&np.all(upper>=point-radius,axis=1))
    edges=np.unique(np.sort(np.vstack([F[selected][:,[0,1]],F[selected][:,[1,2]],F[selected][:,[2,0]]]),axis=1),axis=0)
    crossings=[]
    for edge in edges:
        p,q=V[edge];a=np.minimum(p,q);b=np.maximum(p,q)
        possible=selected[np.all(lower[selected]<=b+1e-9,axis=1)&np.all(upper[selected]>=a-1e-9,axis=1)]
        for index in possible:
            if np.any(np.isin(F[index],edge)):continue
            hit=segment_facet(p,q,tri[index])
            if hit is not None:
                crossings.append(dict(edge_vertex_indices=edge.tolist(),face_index=int(index),
                    edge_points_local_mm=V[edge].tolist(),facet_points_local_mm=tri[index].tolist(),intersection=hit))
    unique={}
    for crossing in crossings:
        key=tuple(np.round(crossing['intersection']['point_local_mm'],8));unique.setdefault(key,crossing)
    used=np.unique(F[selected]);mapping=np.full(len(V),-1,dtype=np.int64);mapping[used]=np.arange(len(used))
    mesh_geom.write_stl(folder/'failure_patch_OPEN.stl',V[used],mapping[F[selected]],header='OPEN diagnostic source surface patch, not printable or a volume')
    return dict(engine_point_local_mm=point.tolist(),source_coordinate_point_mm=(point+data['translation_source_mm']).tolist(),
        localization_radius_mm=radius,selected_surface_faces=len(selected),candidate_unique_edges=len(edges),
        strict_nonadjacent_segment_facet_crossings=list(unique.values()),
        open_patch=(folder/'failure_patch_OPEN.stl').relative_to(A.ROOT).as_posix(),
        interpretation='local floating-point intersection diagnostic; open patch is not a replacement muscle or full self-intersection proof')

def cube_control():
    import gmsh
    solid=A.I._manifold().Manifold.cube([10.,10.,10.]);V,F=mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
    gmsh.initialize(['gmsh','-nopopup'],readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('control_cube')
        _,volume=M.import_boundary(gmsh,V,F);gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1);gmsh.option.setNumber('Mesh.MeshSizeMax',2.)
        gmsh.model.mesh.generate(3)
        tags,coordinates,_=gmsh.model.mesh.getNodes();order=np.argsort(tags);tags=np.asarray(tags)[order];nodes=np.asarray(coordinates).reshape(-1,3)[order]
        types,_,elements=gmsh.model.mesh.getElements(3,volume)
        if list(types)!=[4]:raise ValueError('Control did not produce linear tetrahedra')
        tets=np.searchsorted(tags,np.asarray(elements[0]).reshape(-1,4));metrics,volumes,quality=M.tetra_metrics(nodes,tets)
        if metrics['negative_or_zero_volume_count'] or abs(metrics['volume_mm3']-1000.)>1e-6:raise ValueError('Meshing API control failed')
        np.savez_compressed(BASE/'control_cube.npz',nodes_local_mm=nodes,tetrahedra=tets,signed_volumes_mm3=volumes,mean_ratio=quality)
        return dict(passed=True,expected_volume_mm3=1000.,metrics=metrics,proves_complete_muscle_mesh_quality=False)
    finally:gmsh.finalize()

def run(replay):
    BASE.mkdir(parents=True,exist_ok=True);diagnostics=[]
    for module in ['infraspinatus','subscapularis']:
        folder=BASE/module/'validated_manifold_h0.8_facetangle1e-06_algorithm10'
        if replay:
            for name in ['mesh_attempt.json','gmsh_log.txt']:
                old=folder/name;backup=folder/('before_stdout_capture_'+name)
                if old.exists() and not backup.exists():shutil.copy2(old,backup)
            code="from pathlib import Path; import sys,runpy; sys.path.insert(0,str(Path.cwd()/'.deps')); sys.argv=['mesh','--module',"+repr(module)+",'--mode','validated_manifold','--facet-angle-tol','0.000001','--algorithm','10']; runpy.run_path('mujoco_arm/tools/build_shoulder_tpu_volume_mesh.py',run_name='__main__')"
            with (folder/'engine_stdout.txt').open('w',encoding='utf8') as stream:
                process=subprocess.run([sys.executable,'-c',code],cwd=A.ROOT,stdout=stream,stderr=subprocess.STDOUT)
            attempt=json.loads((folder/'mesh_attempt.json').read_text(encoding='utf8'))
            if process.returncode not in [0,2] or process.returncode==2 and attempt['geometry_tetrahedralization_pass'] and attempt['numerical_elasticity_mesh_quality_pass']:
                raise ValueError('Mesher subprocess exit status disagrees with its gate report')
        raw=(folder/'engine_stdout.txt').read_text(encoding='utf8',errors='replace')
        match=re.search(r'at point \(([-+\d.eE]+),([-+\d.eE]+),([-+\d.eE]+)\)',raw)
        if not match:raise ValueError('No captured native PLC failure point')
        point=np.array([float(v) for v in match.groups()]);result=localize(folder,point)
        result['module']=module;result['native_stdout_sha256']=A.hashed(folder/'engine_stdout.txt')
        routes=json.loads((A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json').read_text(encoding='utf8'))['routes']
        route=next(r for r in routes if r['module']==module);source_point=np.array(result['source_coordinate_point_mm'])
        result['port_context']={}
        for side,port in route['ports'].items():
            delta=source_point-np.array(port['inner_face_source_mm']);normal=np.array(port['normal_source'])
            result['port_context'][side]=dict(distance_to_inner_face_point_mm=float(np.linalg.norm(delta)),
                normal_offset_mm=float(delta@normal),in_plane_radius_mm=float(np.linalg.norm(delta-normal*(delta@normal))))
        normal=np.array(route['ports']['distal']['normal_source']);tangent=np.array(route['distal_terminal_override']['tangent_source'])
        tangent=tangent-normal*(tangent@normal);tangent/=np.linalg.norm(tangent)
        basis=np.column_stack([tangent,np.cross(normal,tangent),normal])
        result['distal_interface_frame_coordinates_mm']=(basis.T@(source_point-np.array(route['ports']['distal']['inner_face_source_mm']))).tolist()
        infill=json.loads((A.DESIGN/'s1_y_muscle_pair/transition_stock_comparison/transition_trials.json').read_text(encoding='utf8'))
        changed=next(r for r in infill['results'] if r['module']==module and r['variant']=='closed_pockets_infilled')
        boxes=[np.array(c['bounds_source_mm']) for c in changed['explicit_infill_changes']]
        result['engine_point_inside_any_Y2_infill_bounds']=any(np.all(source_point>=b[:3]) and np.all(source_point<=b[3:]) for b in boxes)
        diagnostics.append(result)
        print(module,'localized faces',result['selected_surface_faces'],'strict crossings',len(result['strict_nonadjacent_segment_facet_crossings']),flush=True)
    control=cube_control();attempts=[]
    for f in sorted(BASE.rglob('mesh_attempt.json')):
        data=json.loads(f.read_text(encoding='utf8'));attempts.append(dict(file=f.relative_to(A.ROOT).as_posix(),sha256=A.hashed(f),
            module=data['module'],mode=data['mode'],gmsh_algorithm_3d=data.get('gmsh_algorithm_3d',1),
            geometry_pass=data['geometry_tetrahedralization_pass'],elasticity_quality_pass=data['numerical_elasticity_mesh_quality_pass'],error=data.get('error')))
    output=dict(scope=__doc__,script_sha256=A.hashed(Path(__file__).resolve()),mesher_script_sha256=A.hashed(Path(M.__file__).resolve()),
        control=control,localized_failures=diagnostics,attempts=attempts,
        complete_muscle_meshes_ready_for_elasticity=any(r['geometry_pass'] and r['elasticity_quality_pass'] for r in attempts),
        physical_material_calibrated=False,stage_pass=False,actual_drive_range_deg=None)
    (BASE/'diagnostics.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')
    print('Control cube passed; complete TPU elasticity readiness',output['complete_muscle_meshes_ready_for_elasticity'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--replay-hxt',action='store_true')
    args=parser.parse_args();run(args.replay_hxt)
