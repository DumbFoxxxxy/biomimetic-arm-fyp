"""Traceable tetrahedralization attempts for current complete TPU bellies.

Reads the real STL in Python, records exact and 1e-4mm weld topology, and
imports the validated boundary into Gmsh without native-path workarounds.
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

def run(module,size,mode):
    import gmsh
    selection_file=ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=json.loads(selection_file.read_text(encoding='utf8'))
    assembly_file=ROOT/selection['assembly_report'];assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    if A.hashed(assembly_file)!=selection['assembly_sha256']:raise ValueError('Stale candidate')
    source=assembly['sources']['tpu:'+module];solid,_,_=A.load_checked(ROOT/source['file'],source['translation_source_mm'],source['sha256'])
    soup,faces=mesh_geom.read_stl(ROOT/source['file']);exact,inv=np.unique(soup,axis=0,return_inverse=True);raw_faces=inv[faces]
    V,F=mesh_geom.weld(soup,faces,decimals=4)
    _,boundary,nonmanifold=mesh_geom.edge_census(F)
    if boundary or nonmanifold:raise ValueError('Current validated boundary is not closed')
    # Account for the existing geometric import convention. Rounded keys only
    # choose the retained first real coordinate; no vertex grid-snapping occurs.
    _,first,inverse=np.unique(np.round(soup,4),axis=0,return_index=True,return_inverse=True)
    moved=np.linalg.norm(soup-soup[first[inverse]],axis=1)
    dest=assembly_file.parent/'volume_mesh'/module/(mode+'_h'+f'{size:g}');dest.mkdir(parents=True,exist_ok=True)
    raw_census=mesh_geom.edge_census(raw_faces)
    tri=V[F];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2.
    metadata=dict(scope=__doc__,module=module,assembly_source=assembly_file.relative_to(ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        source=source,gmsh_version=gmsh.__version__,script_sha256=A.hashed(Path(__file__).resolve()),
        mode=mode,requested_max_interior_size_mm=size,
        input_preprocessing=dict(raw_unique_vertices=len(exact),raw_triangles=len(raw_faces),raw_boundary_edges=raw_census[1],raw_nonmanifold_edges=raw_census[2],
            validation_weld_key_decimals=4,validation_weld_key_pitch_mm=.0001,maximum_reassigned_vertex_distance_mm=float(moved.max()),
            preprocessed_vertices=len(V),preprocessed_triangles=len(F),dropped_degenerate_triangles=len(raw_faces)-len(F),
            preprocessed_boundary_edges=boundary,preprocessed_nonmanifold_edges=nonmanifold,
            triangle_area_min_mm2=float(areas.min()),triangle_area_p01_mm2=float(np.quantile(areas,.01)),
            interpretation='same rounded-key merge convention as current mesh checks; explicit analysis preprocessing, source STL unchanged'),
        original_material_volume_mm3=float(solid.volume()),geometry_tetrahedralization_pass=False,numerical_elasticity_mesh_quality_pass=False,
        physical_material_model_calibrated=False,nonlinear_solve_completed=False,stage_pass=False,actual_drive_range_deg=None)
    np.savez_compressed(dest/'input_boundary.npz',nodes_local_mm=V,triangles=F,translation_source_mm=np.array(source['translation_source_mm']))
    (dest/'mesh_attempt.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
    gmsh.initialize(['gmsh','-nopopup'],readConfigFiles=False);gmsh.logger.start();started=time.monotonic()
    try:
        gmsh.option.setNumber('General.Terminal',1);gmsh.option.setNumber('General.NumThreads',2)
        gmsh.model.add('current_TPU_'+module)
        surface,volume=import_boundary(gmsh,V,F)
        gmsh.option.setNumber('Mesh.MeshOnlyEmpty',1)
        gmsh.option.setNumber('Mesh.Algorithm3D',1)
        gmsh.option.setNumber('Mesh.MeshSizeMax',size)
        gmsh.option.setNumber('Mesh.MeshSizeMin',min(size,.2))
        gmsh.option.setNumber('Mesh.Optimize',1)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
        gmsh.model.mesh.generate(3)
        tags,coords,_=gmsh.model.mesh.getNodes();order=np.argsort(tags);sorted_tags=np.asarray(tags)[order];nodes=np.asarray(coords).reshape(-1,3)[order]
        types,element_tags,elements=gmsh.model.mesh.getElements(3,volume)
        if list(types)!=[4]:raise ValueError('Expected only linear four-node tetrahedra')
        tetra=np.searchsorted(sorted_tags,np.asarray(elements[0]).reshape(-1,4))
        metrics,volumes,quality=tetra_metrics(nodes,tetra)
        np.savez_compressed(dest/'volume_mesh.npz',nodes_local_mm=nodes,tetrahedra=tetra,translation_source_mm=np.array(source['translation_source_mm']),
            signed_volumes_mm3=volumes,mean_ratio=quality)
        delta=abs(metrics['volume_mm3']-solid.volume());metadata.update(metrics=metrics,volume_difference_mm3=float(delta),
            geometry_tetrahedralization_pass=metrics['negative_or_zero_volume_count']==0 and delta/max(1.,solid.volume())<1e-5,
            numerical_elasticity_mesh_quality_pass=metrics['negative_or_zero_volume_count']==0 and metrics['mean_ratio_below_0_05_count']==0,
            mesh_file=(dest/'volume_mesh.npz').relative_to(ROOT).as_posix(),mesh_sha256=A.hashed(dest/'volume_mesh.npz'))
    except Exception as exc:
        metadata['error_type']=type(exc).__name__;metadata['error']=str(exc)
        print('Mesh attempt rejected:',module,type(exc).__name__,str(exc),flush=True)
    finally:
        metadata['elapsed_s']=time.monotonic()-started
        (dest/'gmsh_log.txt').write_text('\n'.join(gmsh.logger.get())+'\n',encoding='utf8')
        gmsh.logger.stop();gmsh.finalize()
        metadata['source_bytes_unchanged']=A.hashed(ROOT/source['file'])==source['sha256']
        (dest/'mesh_attempt.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
    print(module,'geometry',metadata['geometry_tetrahedralization_pass'],'quality',metadata['numerical_elasticity_mesh_quality_pass'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--module',choices=['infraspinatus','subscapularis'],default='infraspinatus')
    parser.add_argument('--mesh-size',type=float,default=.8);parser.add_argument('--mode',choices=['preserve'],default='preserve')
    args=parser.parse_args();run(args.module,args.mesh_size,args.mode)
