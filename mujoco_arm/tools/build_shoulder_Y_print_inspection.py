"""Crop source-faithful inspection segments; no mechanical or printing pass.

Segments have flat artificial cuts, no grips and no load calibration. They
test fabrication details only, never substitute for whole-muscle tests.
"""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom

DEST = A.DESIGN / 's1_y_muscle_pair/print_inspection'

def write_readback(solid, path):
    v, f = mesh_geom.manifold_mesh_to_arrays(solid.to_mesh())
    mesh_geom.write_stl(str(path), v, f, header='Source-cropped inspection only; no drive or printing pass')
    vertices, faces = mesh_geom.read_stl(str(path))
    vertices, faces = mesh_geom.weld(vertices, faces)
    _, boundary, nonmanifold = mesh_geom.edge_census(faces)
    restored = A.I._manifold().Manifold(A.I._manifold().Mesh64(vertices, faces.astype(np.uint32)))
    components = [float(p.volume()) for p in restored.decompose()]
    return restored, dict(file=path.relative_to(A.ROOT).as_posix(), sha256=A.hashed(path),
        boundary_edges=boundary, nonmanifold_edges=nonmanifold,
        component_volumes_mm3=components, positive_components=sum(v > 1e-6 for v in components),
        relative_volume_roundtrip_error=abs(restored.volume()-solid.volume())/max(1.,solid.volume()),
        no_material_removed_or_added_as_repair=True)

def run():
    assembly_file=A.DESIGN/'s1_y_muscle_pair/assembly_check.json'
    generation_file=assembly_file.parent/'candidate_generation.json'
    assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    generation=json.loads(generation_file.read_text(encoding='utf8'))
    if A.hashed(generation_file)!=assembly['generation_sha256']:
        raise ValueError('Stale generation report')
    DEST.mkdir(parents=True, exist_ok=True)
    specimens=[]
    for row in generation['attempts']:
        if not row.get('selected'): continue
        source=assembly['sources']['tpu:'+row['module']]
        solid,_,_=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])
        frames=row['transported_frames']['stations']
        nodes=row['transported_frames']['nodes']
        # Exclude terminal stations so the two inspection regions cover actual
        # belly details rather than mostly solid mounting interfaces.
        sizes=[]
        for i in range(1,len(frames)-1):
            pts=np.array([p['source_mm'] for p in nodes if p['index'][2]==i])
            F=np.array(frames[i]['paper_basis_columns_source'])
            centre=np.array(frames[i]['centre_source_mm'])
            span=np.ptp((pts-centre)@F,axis=0)
            sizes.append((float(span[0]*span[2]),i,span))
        chosen=[('narrow',min(sizes,key=lambda v:v[0])),('broad',max(sizes,key=lambda v:v[0]))]
        for kind,(area,i,span) in chosen:
            centre=np.array(frames[i]['centre_source_mm']);F=np.array(frames[i]['paper_basis_columns_source'])
            half=float(row['transported_frames']['station_spacing_mm'])*1.05
            affine=np.c_[F.T,-F.T@centre]
            local=solid.transform(affine)
            box=A.I._manifold().Manifold.cube([400.,2*half,400.]).translate([-200.,-half,-200.])
            cropped=local ^ box
            if cropped.is_empty(): raise ValueError('Empty crop')
            # Rotate local Y to print-candidate +Z, keeping a proper rotation.
            R=np.array([[1.,0.,0.],[0.,0.,-1.],[0.,1.,0.]])
            candidate=cropped.transform(np.c_[R,np.zeros(3)])
            bounds=np.array(candidate.bounding_box()).reshape(2,3)
            shift=np.array([-(bounds[0,0]+bounds[1,0])/2.,-(bounds[0,1]+bounds[1,1])/2.,-bounds[0,2]])
            candidate=candidate.translate(shift)
            path=DEST/(row['module']+'_'+kind+'.stl')
            restored,mesh=write_readback(candidate,path)
            v,f=mesh_geom.manifold_mesh_to_arrays(restored.to_mesh64())
            triangles=v[f];cross=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
            lengths=np.linalg.norm(cross,axis=1);areas=lengths/2.
            normals=cross/np.maximum(lengths[:,None],1e-30)
            bottom=np.max(triangles[:,:,2],axis=1)<.01
            # 45 degrees is a geometric inspection bin, not a selected FDM
            # capability. Exclude artificial bottom faces from the area bin.
            down=(normals[:,2]<-np.cos(np.pi/4)) & ~bottom
            total=float(areas.sum())
            T=R@F.T;t=shift-T@centre
            components=mesh['positive_components']
            specimens.append(dict(module=row['module'],region=kind,station=i,
                selection_rule='minimum/maximum mapped X*Z grid span among nonterminal stations',
                nominal_mapped_section_area_proxy_mm2=area,mapped_X_Y_Z_span_mm=span.tolist(),
                source_mesh=source['file'],source_sha256=source['sha256'],
                cut_planes_in_source=dict(centre_mm=centre.tolist(),normal=F[:,1].tolist(),half_separation_mm=half),
                source_to_candidate_affine=np.c_[T,t].tolist(),
                candidate_to_source_affine=np.c_[T.T,-T.T@t].tolist(),
                print_candidate_axis='local station Y to +Z, unselected manufacturing orientation',
                mesh=mesh,dimensions_mm=(np.array(restored.bounding_box())[3:]-np.array(restored.bounding_box())[:3]).tolist(),
                downward_surface_area_beyond_45deg_mm2=float(areas[down].sum()),
                downward_surface_fraction_excluding_bottom=float(areas[down].sum()/max(total,1e-30)),
                overhang_interpretation='triangle-normal diagnostic only; not support, bridging or slicer analysis',
                single_material_component=components==1,
                geometry_readback_usable=(mesh['boundary_edges']==mesh['nonmanifold_edges']==0 and components==1 and mesh['relative_volume_roundtrip_error']<1e-5),
                print_process_verified=False,channel_threading_verified=False,mechanical_calibration_usable=False,
                exclusions=['artificial section cuts replace whole-muscle boundary conditions','no grips or load transfer interfaces','cropped channel has artificial openings','no actual slicing or material test']))
            print(row['module'],kind,'components',components,'boundary',mesh['boundary_edges'],'nonmanifold',mesh['nonmanifold_edges'],flush=True)
    report=dict(scope=__doc__,assembly_source=assembly_file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        generation_sha256=A.hashed(generation_file),script_sha256=A.hashed(Path(__file__).resolve()),
        source_meshes_unchanged=True,specimens=specimens,actual_drive_range_deg=None,stage_pass=False)
    (DEST/'inspection_segments.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

if __name__=='__main__': run()
