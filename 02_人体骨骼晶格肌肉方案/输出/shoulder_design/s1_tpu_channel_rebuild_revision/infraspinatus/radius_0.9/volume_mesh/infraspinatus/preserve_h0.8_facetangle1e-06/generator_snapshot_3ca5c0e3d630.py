"""Traceable tetrahedralization attempts for current complete TPU bellies.

Reads source STL or explicitly indexed float64 OBJ in Python, records exact
and declared-import topology, and imports the boundary into Gmsh without
native-path workarounds. Legacy STL retains its existing 1e-4mm key weld.
Volume meshes are analysis inputs only, never calibrated TPU or drive proof.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom

ROOT=A.ROOT
VENDOR=ROOT/'.deps/gmsh_vendor'
if VENDOR.exists():sys.path.insert(0,str(VENDOR))

def import_boundary(gmsh,V,F):
    surface=gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2,surface,np.arange(1,len(V)+1,dtype=np.uint64),V.ravel())
    gmsh.model.mesh.addElementsByType(surface,2,np.arange(1,len(F)+1,dtype=np.uint64),(F+1).astype(np.uint64).ravel())
    loop=gmsh.model.geo.addSurfaceLoop([surface]);volume=gmsh.model.geo.addVolume([loop]);gmsh.model.geo.synchronize()
    return surface,volume

def tetra_metrics(nodes,tets):
    p=nodes[tets];det=np.einsum('ij,ij->i',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))
    volumes=det/6.;edge2=sum(np.sum((p[:,i]-p[:,j])**2,axis=1) for i in range(4) for j in range(i))
    quality=12.*np.cbrt(3.*np.abs(volumes))**2/np.maximum(edge2,1e-300)
    return dict(tetra_count=len(tets),negative_or_zero_volume_count=int(np.sum(volumes<=0)),
        volume_mm3=float(volumes.sum()),minimum_volume_mm3=float(volumes.min()),
        mean_ratio_min=float(quality.min()),mean_ratio_quantiles=np.quantile(quality,[.01,.1,.5]).tolist(),
        mean_ratio_below_0_05_count=int(np.sum(quality<.05))),volumes,quality

def load_boundary(source):
    """Respect the declared import convention; never bucket-weld indexed OBJ."""
    path=ROOT/source['file']
    if A.hashed(path)!=source['sha256']:raise ValueError('Changed mesh input')
    indexed=source.get('import_convention')=='indexed_float64_obj_no_bucket_weld'
    if indexed:
        if path.suffix.lower()!='.obj':raise ValueError('Indexed float64 convention requires OBJ')
        soup,faces=mesh_geom.load_obj(str(path));V,F=soup.copy(),faces.copy()
        _,boundary,nonmanifold=mesh_geom.edge_census(F)
        if boundary or nonmanifold:raise ValueError('Indexed OBJ boundary is open or nonmanifold')
        solid=A.I._manifold().Manifold(A.I._manifold().Mesh64(
            np.ascontiguousarray(V+np.array(source['translation_source_mm']),dtype=np.float64),
            np.ascontiguousarray(F,dtype=np.uint64)))
        if solid.volume()<=0 or str(solid.status()).split('.')[-1]!='NoError':raise ValueError('Invalid indexed solid')
        decimals=pitch=None;max_move=0.
        interpretation='explicit indexed float64 OBJ: no welding, snapping or triangle deletion; exact-coordinate census is diagnostic only'
    else:
        solid,_,_=A.load_checked(path,source['translation_source_mm'],source['sha256'])
        if path.suffix.lower()!='.stl':raise ValueError('Legacy mesh convention expects STL')
        soup,faces=mesh_geom.read_stl(path);V,F=mesh_geom.weld(soup,faces,decimals=4)
        _,boundary,nonmanifold=mesh_geom.edge_census(F)
        _,first,inverse=np.unique(np.round(soup,4),axis=0,return_index=True,return_inverse=True)
        moved=np.linalg.norm(soup-soup[first[inverse]],axis=1)
        decimals=4;pitch=.0001;max_move=float(moved.max())
        interpretation='same rounded-key merge convention as current mesh checks; explicit analysis preprocessing, source STL unchanged'
    exact,inv=np.unique(soup,axis=0,return_inverse=True);raw_faces=inv[faces]
    if boundary or nonmanifold:raise ValueError('Current boundary is not closed')
    raw_census=mesh_geom.edge_census(raw_faces)
    tri=V[F];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2.
    preprocessing=dict(raw_unique_vertices=len(exact),raw_triangles=len(raw_faces),
        raw_boundary_edges=raw_census[1],raw_nonmanifold_edges=raw_census[2],
        validation_weld_key_decimals=decimals,validation_weld_key_pitch_mm=pitch,
        maximum_reassigned_vertex_distance_mm=max_move,preprocessed_vertices=len(V),preprocessed_triangles=len(F),
        dropped_degenerate_triangles=len(raw_faces)-len(F),preprocessed_boundary_edges=boundary,
        preprocessed_nonmanifold_edges=nonmanifold,triangle_area_min_mm2=float(areas.min()),
        triangle_area_p01_mm2=float(np.quantile(areas,.01)),interpretation=interpretation)
    return solid,V,F,preprocessing

def run(module,size,mode,tolerance,facet_angle,algorithm,assembly_override=None):
    import gmsh
    selection_file=ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=json.loads(selection_file.read_text(encoding='utf8'))
    assembly_file=ROOT/selection['assembly_report'] if assembly_override is None else Path(assembly_override).resolve()
    assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    if assembly_override is None and A.hashed(assembly_file)!=selection['assembly_sha256']:raise ValueError('Stale candidate')
    source=assembly['sources']['tpu:'+module];solid,V,F,preprocessing=load_boundary(source)
    suffix=mode+('_e'+f'{tolerance:g}' if mode=='analysis_simplify' else '')+'_h'+f'{size:g}'
    if facet_angle!=.1:suffix+='_facetangle'+f'{facet_angle:g}'
    if algorithm!=1:suffix+='_algorithm'+str(algorithm)
    dest=assembly_file.parent/'volume_mesh'/module/suffix;dest.mkdir(parents=True,exist_ok=True)
    metadata=dict(scope=__doc__,module=module,assembly_source=assembly_file.relative_to(ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        source=source,gmsh_version=gmsh.__version__,script_sha256=A.hashed(Path(__file__).resolve()),explicit_nondefault_assembly=assembly_override is not None,
        mode=mode,requested_max_interior_size_mm=size,gmsh_near_facet_angle_tolerance_deg=facet_angle,gmsh_algorithm_3d=algorithm,
        near_facet_option_is_numerical_recovery_not_an_elasticity_acceptance_gate=True,
        input_preprocessing=preprocessing,
        original_material_volume_mm3=float(solid.volume()),geometry_tetrahedralization_pass=False,numerical_elasticity_mesh_quality_pass=False,
        physical_material_model_calibrated=False,nonlinear_solve_completed=False,stage_pass=False,actual_drive_range_deg=None)
    if mode in ['analysis_simplify','validated_manifold']:
        origin=np.array(source['translation_source_mm'])
        approximate=solid.translate(-origin)
        if mode=='analysis_simplify':approximate=approximate.simplify(tolerance)
        V,F=mesh_geom.manifold_mesh_to_arrays(approximate.to_mesh64())
        restored=approximate.translate(origin)
        _,bnd,nm=mesh_geom.edge_census(F)
        metadata['analysis_surface']=dict(simplification_tolerance_mm=tolerance if mode=='analysis_simplify' else None,vertices=len(V),triangles=len(F),
            boundary_edges=bnd,nonmanifold_edges=nm,
            added_material_signed_mm3=float((restored-solid).volume()),removed_material_signed_mm3=float((solid-restored).volume()),
            material_volume_mm3=float(restored.volume()),
            source_geometry_replaced=False,geometry_approximation_not_a_collision_tolerance=True)
        if bnd or nm:raise ValueError('Analysis simplification failed topology preflight')
    metadata['gmsh_python_sha256']=A.hashed(Path(gmsh.__file__).resolve())
    if (VENDOR/'gmsh-4.15.dll').exists():metadata['gmsh_windows_dll_sha256']=A.hashed(VENDOR/'gmsh-4.15.dll')
    np.savez_compressed(dest/'input_boundary.npz',nodes_local_mm=V,triangles=F,translation_source_mm=np.array(source['translation_source_mm']))
    (dest/'mesh_attempt.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
    gmsh.initialize(['gmsh','-nopopup'],readConfigFiles=False);gmsh.logger.start();started=time.monotonic()
    try:
        gmsh.option.setNumber('General.Terminal',1);gmsh.option.setNumber('General.NumThreads',2)
        gmsh.model.add('current_TPU_'+module)
        surface,volume=import_boundary(gmsh,V,F)
        gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)
        gmsh.option.setNumber('Mesh.Algorithm3D',algorithm)
        gmsh.option.setNumber('Mesh.MeshSizeMax',size)
        gmsh.option.setNumber('Mesh.MeshSizeMin',min(size,.2))
        gmsh.option.setNumber('Mesh.Optimize',1)
        gmsh.option.setNumber('Mesh.AngleToleranceFacetOverlap',facet_angle)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
        gmsh.model.mesh.generate(3)
        tags,coords,_=gmsh.model.mesh.getNodes();order=np.argsort(tags);sorted_tags=np.asarray(tags)[order];nodes=np.asarray(coords).reshape(-1,3)[order]
        types,element_tags,elements=gmsh.model.mesh.getElements(3,volume)
        if list(types)!=[4]:raise ValueError('Expected only linear four-node tetrahedra')
        raw_connectivity=np.asarray(elements[0]).reshape(-1,4)
        np.savez_compressed(dest/'raw_engine_connectivity.npz',node_tags=sorted_tags,nodes_local_mm=nodes,
            tetrahedron_node_tags=raw_connectivity,element_tags=np.asarray(element_tags[0]))
        if len(np.unique(sorted_tags))!=len(sorted_tags):raise ValueError('Engine returned duplicate node tags')
        tetra=np.searchsorted(sorted_tags,raw_connectivity)
        if np.any(tetra>=len(sorted_tags)) or not np.array_equal(sorted_tags[tetra],raw_connectivity):
            raise ValueError('Element references absent node tags; refusing approximate tag lookup')
        metadata['raw_engine_connectivity_file']=(dest/'raw_engine_connectivity.npz').relative_to(ROOT).as_posix()
        metadata['raw_engine_connectivity_sha256']=A.hashed(dest/'raw_engine_connectivity.npz')
        metadata['all_element_node_tags_exactly_resolved']=True
        metrics,volumes,quality=tetra_metrics(nodes,tetra)
        np.savez_compressed(dest/'volume_mesh.npz',nodes_local_mm=nodes,tetrahedra=tetra,translation_source_mm=np.array(source['translation_source_mm']),
            signed_volumes_mm3=volumes,mean_ratio=quality)
        target_volume=metadata.get('analysis_surface',{}).get('material_volume_mm3',float(solid.volume()))
        delta=abs(metrics['volume_mm3']-target_volume);metadata.update(metrics=metrics,volume_difference_from_analysis_boundary_mm3=float(delta),
            volume_difference_from_original_mm3=float(metrics['volume_mm3']-solid.volume()),
            geometry_tetrahedralization_pass=metrics['negative_or_zero_volume_count']==0 and delta/max(1.,target_volume)<1e-5,
            numerical_elasticity_mesh_quality_pass=metrics['negative_or_zero_volume_count']==0 and metrics['mean_ratio_below_0_05_count']==0,
            mesh_file=(dest/'volume_mesh.npz').relative_to(ROOT).as_posix(),mesh_sha256=A.hashed(dest/'volume_mesh.npz'))
    except Exception as exc:
        metadata['error_type']=type(exc).__name__;metadata['error']=str(exc)
        try:
            tags,coordinates,_=gmsh.model.mesh.getNodes()
            bad=gmsh.model.mesh.getLastNodeError()
            metadata['reported_bad_node_tags']=list(map(int,bad))
            metadata['reported_bad_entities']=[list(map(int,e)) for e in gmsh.model.mesh.getLastEntityError()]
            np.savez_compressed(dest/'failed_mesher_nodes.npz',node_tags=np.asarray(tags),nodes_local_mm=np.asarray(coordinates).reshape(-1,3))
        except Exception as diagnostic_error:metadata['diagnostic_capture_error']=str(diagnostic_error)
        print('Mesh attempt rejected:',module,type(exc).__name__,str(exc),flush=True)
    finally:
        metadata['elapsed_s']=time.monotonic()-started
        (dest/'gmsh_log.txt').write_text('\n'.join(gmsh.logger.get())+'\n',encoding='utf8')
        gmsh.logger.stop();gmsh.finalize()
        metadata['source_bytes_unchanged']=A.hashed(ROOT/source['file'])==source['sha256']
        (dest/'mesh_attempt.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
    print(module,'geometry',metadata['geometry_tetrahedralization_pass'],'quality',metadata['numerical_elasticity_mesh_quality_pass'],flush=True)
    return metadata

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--module',choices=['infraspinatus','subscapularis'],default='infraspinatus')
    parser.add_argument('--mesh-size',type=float,default=.8);parser.add_argument('--mode',choices=['preserve','validated_manifold','analysis_simplify'],default='preserve')
    parser.add_argument('--simplify-tolerance',type=float,default=.005)
    parser.add_argument('--facet-angle-tol',type=float,default=.1)
    parser.add_argument('--algorithm',type=int,choices=[1,10],default=1)
    parser.add_argument('--assembly',type=Path,help='Explicit research candidate; leaves the default selection unchanged')
    args=parser.parse_args();result=run(args.module,args.mesh_size,args.mode,args.simplify_tolerance,args.facet_angle_tol,args.algorithm,args.assembly)
    if not result['geometry_tetrahedralization_pass'] or not result['numerical_elasticity_mesh_quality_pass']:raise SystemExit(2)
