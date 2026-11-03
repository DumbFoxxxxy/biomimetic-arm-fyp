"""Independently replay source crops and check unchanged current-muscle meshes."""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom

DEST=A.DESIGN/'s1_y_muscle_pair/print_inspection'

def run():
    source=DEST/'inspection_segments.json'
    report=json.loads(source.read_text(encoding='utf8'))
    assembly=json.loads((A.ROOT/report['assembly_source']).read_text(encoding='utf8'))
    if A.hashed(A.ROOT/report['assembly_source'])!=report['assembly_sha256']:raise ValueError('Stale assembly')
    if A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_Y_print_inspection.py')!=report['script_sha256']:raise ValueError('Changed generator')
    results=[];originals={};errors=[];void_audits=[]
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes']
    for row in report['specimens']:
        entry=assembly['sources']['tpu:'+row['module']]
        if row['source_sha256']!=entry['sha256']:raise ValueError('Mismatched muscle')
        if row['module'] not in originals:
            full,_,_=A.load_checked(A.ROOT/entry['file'],entry['translation_source_mm'],entry['sha256'])
            originals[row['module']]=full
            head={'infraspinatus':'FJ1500','subscapularis':'FJ1504'}[row['module']]
            route=next(r for r in routes if r['head']==head)
            wire=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3)
            lumen=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.86)
            for shell in full.decompose():
                if shell.volume()>=-1e-6: continue
                v,f=mesh_geom.manifold_mesh_to_arrays(shell.to_mesh64())
                # Reverse a copy solely to measure the EMPTY cavity volume.
                # Never union this positive probe into the source muscle.
                probe=A.I._manifold().Manifold(A.I._manifold().Mesh64(
                    np.array(v,dtype=np.float64,order='C',copy=True),
                    np.array(f[:,::-1],dtype=np.uint64,order='C',copy=True)))
                if probe.volume()<=0:raise ValueError('Invalid positive cavity probe')
                void_audits.append(dict(module=row['module'],empty_volume_mm3=-float(shell.volume()),
                    bounds_source_mm=list(shell.bounding_box()),
                    overlaps_nominal_wire_mm3=A.I.intersection_volume(probe,wire),
                    overlaps_nominal_routed_lumen_mm3=A.I.intersection_volume(probe,lumen),
                    cavity_filled_in_design=False,physical_printability_verified=False))
        full=originals[row['module']]
        mesh=row['mesh'];actual,_,meta=A.load_checked(A.ROOT/mesh['file'],[0.,0.,0.],mesh['sha256'])
        inverse=np.array(row['candidate_to_source_affine']);T=np.array(row['source_to_candidate_affine'])
        source_crop=actual.transform(inverse)
        # Reconstruct cut slab from the recorded source planes, independently
        # of the generator's crop-local intermediate geometry.
        candidate_full=full.transform(T)
        c=np.array(row['cut_planes_in_source']['centre_mm']);n=np.array(row['cut_planes_in_source']['normal'])
        point=T[:,:3]@c+T[:,3];half=row['cut_planes_in_source']['half_separation_mm']
        low=point[2]-half
        expected=candidate_full ^ A.I._manifold().Manifold.cube([800.,800.,2*half]).translate([-400.,-400.,low])
        extra=float((actual-expected).volume());lost=float((expected-actual).volume())
        allowance=max(.001,float(expected.volume())*1e-5)
        cavities=[float(p.volume()) for p in actual.decompose() if p.volume() < -1e-6]
        original_cavities=[float(p.volume()) for p in full.decompose() if p.volume() < -1e-6]
        out=dict(module=row['module'],region=row['region'],source_sha256=entry['sha256'],
            readback_boundary_edges=meta['boundary_edges'],readback_nonmanifold_edges=meta['nonmanifold_edges'],
            replay_extra_material_mm3=extra,replay_missing_material_mm3=lost,report_volume_limit_mm3=allowance,
            inverse_rotation_orthogonality_error=float(np.max(abs(inverse[:,:3].T@inverse[:,:3]-np.eye(3)))),
            candidate_positive_component_volumes_mm3=[float(p.volume()) for p in actual.decompose() if p.volume()>1e-6],
            candidate_enclosed_void_shell_volumes_mm3=cavities,
            whole_source_enclosed_void_shell_volumes_mm3=original_cavities,
            crop_replay_consistent=extra<allowance and lost<allowance,
            printing_or_drive_verified=False)
        if not out['crop_replay_consistent']:errors.append(row['module']+':'+row['region'])
        results.append(out)
        print(row['module'],row['region'],'extra/lost',extra,lost,'source_voids',original_cavities,flush=True)
    output=dict(scope=__doc__,source_sha256=A.hashed(source),assembly_sha256=report['assembly_sha256'],
        script_sha256=A.hashed(Path(__file__).resolve()),routes_source=route_file.relative_to(A.ROOT).as_posix(),routes_sha256=A.hashed(route_file),
        specimens=results,whole_source_void_audits=void_audits,errors=errors,
        stage_pass=False,actual_drive_range_deg=None)
    (DEST/'readback_checks.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    if errors:raise ValueError('Crop replay mismatch')

if __name__=='__main__':run()
