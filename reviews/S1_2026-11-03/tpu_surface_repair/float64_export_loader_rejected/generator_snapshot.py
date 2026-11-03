"""Bounded, additive junction repair trials driven by whole-surface crossings.

Complete original muscles remain references. Local radius-0.05mm stock is
added at independently verified crossings, new stock is kept out of nominal lumens, and
every material change/export is reported. No default selection changes.
"""
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from diagnose_shoulder_tpu_volume_mesh import segment_facet
from build_shoulder_Y_print_inspection import write_readback
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_tpu_surface_global_revision'

def scan(solid):
    V,F=mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
    # Reset source face-property IDs for unambiguous native BVH lookup after
    # Boolean unions of many independently generated patch bodies.
    prepared=A.I._manifold().Manifold(A.I._manifold().Mesh64(
        np.array(V,dtype=np.float64,order='C',copy=True),np.array(F,dtype=np.uint64,order='C',copy=True)))
    mesh=prepared.to_mesh64();V,F=mesh_geom.manifold_mesh_to_arrays(mesh);ids=np.array(mesh.face_id,dtype=np.int64);lookup={}
    for i,face in enumerate(ids):lookup.setdefault(int(face),[]).append(i)
    edges=np.unique(np.sort(np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0);found={}
    for edge in edges:
        p,q=V[edge];seen=set()
        for hit in prepared.ray_cast(p,q):
            for i in lookup.get(int(hit.face_id),[]):
                if i in seen:continue
                seen.add(i)
                if np.any(F[i]==edge[0]) or np.any(F[i]==edge[1]):continue
                crossing=segment_facet(p,q,V[F[i]])
                if crossing is None:continue
                point=np.array(crossing['point_local_mm']);key=tuple(np.round(point,6))
                if key not in found:
                    found[key]=dict(point_source_mm=point.tolist(),edge_source_mm=V[edge].tolist(),facet_source_mm=V[F[i]].tolist(),
                        parameters=[crossing['t'],crossing['u'],crossing['v']])
    return dict(edges_scanned=len(edges),triangles=len(F),count=len(found),crossings=list(found.values()),
        face_id_reset_volume_difference_mm3=float(prepared.volume()-solid.volume()))

def write_obj(solid,path,origin):
    V,F=mesh_geom.manifold_mesh_to_arrays(solid.translate(-origin).to_mesh64())
    with path.open('w',encoding='utf8') as stream:
        stream.write('# Complete TPU research repair, explicit additive stock; no drive or manufacturing pass\n')
        for v in V:stream.write('v %.17g %.17g %.17g\n'%tuple(v))
        for f in F+1:stream.write('f %d %d %d\n'%tuple(f))
    actual,_,meta=A.load_checked(path,origin,A.hashed(path));meta['translation_source_mm']=origin.tolist()
    return actual,meta

def run():
    selection=json.loads((A.ROOT/'design/S1_CURRENT_CANDIDATE.json').read_text(encoding='utf8'))
    parent_file=A.ROOT/selection['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=selection['assembly_sha256']:raise ValueError('Changed source selection')
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in parent['sources'].items()}
    routes=json.loads((A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json').read_text(encoding='utf8'))['routes']
    DEST.mkdir(parents=True,exist_ok=True);results=[]
    for module in ['infraspinatus','subscapularis']:
        old=S['tpu:'+module];current=old;source=parent['sources']['tpu:'+module];origin=np.array(source['translation_source_mm'])
        r=next(r for r in routes if r['module']==module);lumen=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.86)
        rounds=[];folder=DEST/module;folder.mkdir(exist_ok=True)
        for iteration in range(5):
            audit=scan(current);entry=dict(iteration=iteration,surface_audit=audit)
            rounds.append(entry)
            print(module,'global round',iteration,'crossings',audit['count'],flush=True)
            if audit['count']==0 or iteration==4:break
            points=[]
            for c in audit['crossings']:
                point=np.array(c['point_source_mm'])
                if not any(np.linalg.norm(point-p)<.001 for p in points):points.append(point)
            stock=A.I._manifold().Manifold.batch_boolean([A.I._manifold().Manifold.sphere(.05,24).translate(p) for p in points],A.I._manifold().OpType.Add)
            # Preserve the entire original wall. Re-drilling the original
            # float32 wall removed 0.014-0.020mm3 and was rejected. Subtract
            # the nominal lumen from NEW material only, not from the muscle.
            trial=current+(stock-lumen)
            removed=float((old-trial).volume());positive=sum(p.volume()>1e-6 for p in trial.decompose())
            entry['patch']=dict(radius_mm=.05,centres_source_mm=[p.tolist() for p in points],
                added_material_from_parent_signed_mm3=float((trial-old).volume()),removed_parent_material_signed_mm3=removed,
                positive_components=positive)
            if positive!=1 or removed>.001:
                entry['repair_rejected']='material disconnected or original-stock loss exceeds 0.001mm3';break
            current=trial
        exports=[]
        for kind in ['obj64','stl32']:
            if kind=='obj64':actual,meta=write_obj(current,folder/'TPU_global_repair.obj',origin)
            else:
                local,meta=write_readback(current.translate(-origin),folder/'TPU_global_repair.stl');actual=local.translate(origin);meta['translation_source_mm']=origin.tolist()
            after=scan(actual);hits=[];soft=[]
            for key,body in S.items():
                if key=='tpu:'+module:continue
                value=A.I.intersection_volume(actual,body)
                if value>.001:(soft if is_compliant(key) else hits).append(dict(part=key,volume_mm3=value))
            wire_hits=[]
            for route in routes:
                wire=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3);value=A.I.intersection_volume(actual,wire)
                if value>.001:wire_hits.append(dict(head=route['head'],volume_mm3=value))
            exported=dict(kind=kind,mesh=meta,surface_audit=after,rigid_hits_over_threshold=hits,undeformed_TPU_hits_over_threshold=soft,
                active_wire_hits_over_threshold=wire_hits,added_material_signed_mm3=float((actual-old).volume()),
                removed_parent_material_signed_mm3=float((old-actual).volume()),
                nominal_lumen_material_overlap_mm3=A.I.intersection_volume(actual,lumen),adopted_as_current=False,stage_pass=False)
            sources=dict(parent['sources']);sources['tpu:'+module]=meta
            assembly=dict(scope='Complete-context single-muscle global surface repair trial; not a manufacturing or drive pass.',sources=sources,
                parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),source_change=exported,stage_pass=False,actual_drive_range_deg=None)
            (folder/('candidate_assembly_'+kind+'.json')).write_text(json.dumps(assembly,indent=2)+'\n',encoding='utf8')
            exports.append(exported)
            print(module,kind,'readback crossings',after['count'],'rigid hits',len(hits),flush=True)
        row=dict(module=module,rounds=rounds,exports=exports)
        results.append(row);(folder/'repair_history.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf8')
    report=dict(scope=__doc__,parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
        script_sha256=A.hashed(Path(__file__).resolve()),intersection_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/diagnose_shoulder_tpu_volume_mesh.py'),
        export_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_Y_print_inspection.py'),results=results,
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'global_repair_trials.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
