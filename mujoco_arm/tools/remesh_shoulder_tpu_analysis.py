"""Research surface reparametrization following bundled Gmsh tutorial t13.

Replaces the analysis triangulation, not the selected manufacturing surface.
Surface approximation and complete component/contact validation remain gates.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_shoulder_tpu_volume_mesh as B
import mesh_geom
A=B.A


def run(module,assembly_path,size,angle):
    import gmsh
    if assembly_path is None:
        choice=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
        assembly_path=A.ROOT/choice['assembly_report']
        if A.hashed(assembly_path)!=choice['assembly_sha256']:raise ValueError('Changed current assembly')
    assembly_path=assembly_path.resolve();assembly=json.loads(assembly_path.read_text(encoding='utf8'))
    source=assembly['sources']['tpu:'+module];solid,V,F,preprocessing=B.load_boundary(source)
    folder=assembly_path.parent/'surface_remesh'/module/('angle_'+f'{angle:g}'+'_h'+f'{size:g}');folder.mkdir(parents=True,exist_ok=True)
    reference=B.VENDOR/'share/doc/gmsh/tutorials/python/t13.py'
    record=dict(scope=__doc__,module=module,assembly_source=assembly_path.relative_to(A.ROOT).as_posix(),
        assembly_sha256=A.hashed(assembly_path),source=source,script_sha256=A.hashed(Path(__file__).resolve()),
        boundary_loader_sha256=A.hashed(Path(B.__file__)),reference_tutorial=reference.relative_to(A.ROOT).as_posix(),
        reference_tutorial_sha256=A.hashed(reference),input_preprocessing=preprocessing,
        gmsh_version=gmsh.__version__,classification_angle_deg=angle,mesh_size_mm=size,
        geometry_tetrahedralization_pass=False,numerical_elasticity_mesh_quality_pass=False,
        analysis_surface_approximation_verified=False,manufacturing_geometry_replaced=False,
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    started=time.monotonic();gmsh.initialize(['gmsh','-nopopup'],readConfigFiles=False);gmsh.logger.start()
    try:
        gmsh.option.setNumber('General.Terminal',1);gmsh.option.setNumber('General.NumThreads',2)
        gmsh.model.add('TPU_analysis_remeshing_'+module)
        surface=gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2,surface,np.arange(1,len(V)+1,dtype=np.uint64),V.ravel())
        gmsh.model.mesh.addElementsByType(surface,2,np.arange(1,len(F)+1,dtype=np.uint64),(F+1).astype(np.uint64).ravel())
        gmsh.model.mesh.classifySurfaces(math.radians(angle),True,True,math.pi)
        record['classified_entity_counts']={str(dim):len(gmsh.model.getEntities(dim)) for dim in [0,1,2]}
        gmsh.model.mesh.createGeometry()
        surfaces=[tag for _,tag in gmsh.model.getEntities(2)]
        loop=gmsh.model.geo.addSurfaceLoop(surfaces);volume=gmsh.model.geo.addVolume([loop]);gmsh.model.geo.synchronize()
        gmsh.option.setNumber('Mesh.MeshSizeMax',size);gmsh.option.setNumber('Mesh.MeshSizeMin',min(.2,size))
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0);gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0)
        gmsh.option.setNumber('Mesh.Algorithm',6);gmsh.option.setNumber('Mesh.Algorithm3D',10)
        gmsh.option.setNumber('Mesh.Optimize',1)
        gmsh.model.mesh.generate(2)
        tags,xyz,_=gmsh.model.mesh.getNodes();order=np.argsort(tags);tags=np.asarray(tags)[order];nodes=np.asarray(xyz).reshape(-1,3)[order]
        triangles=[]
        for surface in surfaces:
            kind,_,connectivity=gmsh.model.mesh.getElements(2,surface)
            for typ,vertices in zip(kind,connectivity):
                if typ!=2:raise ValueError('Expected only linear boundary triangles')
                triangles.append(np.searchsorted(tags,np.array(vertices).reshape(-1,3)))
        triangles=np.vstack(triangles)
        np.savez_compressed(folder/'remeshed_boundary.npz',nodes_local_mm=nodes,triangles=triangles,
            translation_source_mm=np.array(source['translation_source_mm']))
        _,bnd,nm=mesh_geom.edge_census(triangles)
        record['surface_mesh']=dict(vertices=len(nodes),triangles=len(triangles),boundary_edges=bnd,nonmanifold_edges=nm,
            signed_volume_mm3=mesh_geom.signed_volume(nodes,triangles),
            file=(folder/'remeshed_boundary.npz').relative_to(A.ROOT).as_posix(),sha256=A.hashed(folder/'remeshed_boundary.npz'))
        if bnd or nm:raise ValueError('Remeshed boundary is not closed')
        gmsh.model.mesh.generate(3)
        tags,xyz,_=gmsh.model.mesh.getNodes();order=np.argsort(tags);tags=np.asarray(tags)[order];nodes=np.asarray(xyz).reshape(-1,3)[order]
        kinds,_,connectivity=gmsh.model.mesh.getElements(3,volume)
        if list(kinds)!=[4]:raise ValueError('Expected linear tetrahedra')
        tetra=np.searchsorted(tags,np.asarray(connectivity[0]).reshape(-1,4));metrics,volumes,quality=B.tetra_metrics(nodes,tetra)
        np.savez_compressed(folder/'volume_mesh.npz',nodes_local_mm=nodes,tetrahedra=tetra,
            translation_source_mm=np.array(source['translation_source_mm']),signed_volumes_mm3=volumes,mean_ratio=quality)
        target=abs(record['surface_mesh']['signed_volume_mm3']);delta=abs(metrics['volume_mm3']-target)
        record.update(metrics=metrics,volume_difference_from_remeshed_boundary_mm3=delta,
            volume_difference_from_original_mm3=float(metrics['volume_mm3']-solid.volume()),
            geometry_tetrahedralization_pass=metrics['negative_or_zero_volume_count']==0 and delta/max(1.,target)<1e-5,
            numerical_elasticity_mesh_quality_pass=metrics['negative_or_zero_volume_count']==0 and metrics['mean_ratio_below_0_05_count']==0,
            mesh_file=(folder/'volume_mesh.npz').relative_to(A.ROOT).as_posix(),mesh_sha256=A.hashed(folder/'volume_mesh.npz'))
    except Exception as error:
        record.update(error_type=type(error).__name__,error=str(error))
        print('Remesh rejected',module,type(error).__name__,str(error),flush=True)
    finally:
        record['elapsed_s']=time.monotonic()-started
        (folder/'gmsh_log.txt').write_text('\n'.join(gmsh.logger.get())+'\n',encoding='utf8')
        gmsh.logger.stop();gmsh.finalize()
        record['source_bytes_unchanged']=A.hashed(A.ROOT/source['file'])==source['sha256']
        (folder/'mesh_attempt.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(module,'geometry',record['geometry_tetrahedralization_pass'],'quality',record['numerical_elasticity_mesh_quality_pass'],flush=True)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--module',choices=['infraspinatus','subscapularis'],required=True)
    p.add_argument('--assembly',type=Path);p.add_argument('--mesh-size',type=float,default=.6);p.add_argument('--angle',type=float,default=60.)
    args=p.parse_args();r=run(args.module,args.assembly,args.mesh_size,args.angle)
    if not r['geometry_tetrahedralization_pass'] or not r['numerical_elasticity_mesh_quality_pass']:raise SystemExit(2)
