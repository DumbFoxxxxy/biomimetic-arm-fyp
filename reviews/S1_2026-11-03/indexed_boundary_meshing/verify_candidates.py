"""Independent readback, collapse replay, interface and meshing evidence audit."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'.deps'));sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from build_shoulder_shaft_collar import rod
from check_shoulder_front_liner_keepers import is_compliant
import build_shoulder_tpu_volume_mesh as B

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf8'))
errors=[]
selection_file=ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=read(selection_file)
assembly_file=ROOT/selection['assembly_report'];parent=read(assembly_file)
if sha(assembly_file)!=selection['assembly_sha256']:raise ValueError('Changed selected assembly')
before=read(ROOT/'reviews/S1_2026-11-03/document_optimization/before_SHA256.json')
if sha(selection_file)!=before['design/S1_CURRENT_CANDIDATE.json']:errors.append('Default selection changed')
S={k:A.load_checked(ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in parent['sources'].items()}
routes=read(A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json')['routes']


def raw_indexed(row):
    p=ROOT/row['file']
    if sha(p)!=row['sha256']:raise ValueError('Changed indexed file')
    V,F=mesh_geom.load_obj(str(p));V+=np.array(row['translation_source_mm'])
    _,bnd,nm=mesh_geom.edge_census(F)
    if bnd or nm:errors.append('Bad indexed topology: '+row['file'])
    solid=A.I._manifold().Manifold(A.I._manifold().Mesh64(np.ascontiguousarray(V),np.ascontiguousarray(F,dtype=np.uint64)))
    return V,F,solid


def replay_collapse(original,history,output_V,output_F):
    initial,F=mesh_geom.manifold_mesh_to_arrays(original.to_mesh64());V=initial.copy();F=F.copy()
    members={i:[i] for i in range(len(V))};worst=0.;reversed_faces=0
    for operation in history['operations']:
        a=operation['retained_vertex'];b=operation['removed_vertex'];position=np.array(operation['new_position_source_mm'])
        group=members[a]+members.pop(b);worst=max(worst,float(np.linalg.norm(initial[group]-position,axis=1).max()))
        if len(group)!=operation['original_cluster_vertex_count']:errors.append('Cluster membership mismatch')
        affected=np.flatnonzero(np.any((F==a)|(F==b),axis=1))
        old_vertices=V[F[affected]];changed=F[affected].copy();changed[changed==b]=a
        survives=np.array([len(set(f))==3 for f in changed]);new_vertices=V[changed].copy();new_vertices[changed==a]=position
        old_norm=np.cross(old_vertices[:,1]-old_vertices[:,0],old_vertices[:,2]-old_vertices[:,0])
        new_norm=np.cross(new_vertices[:,1]-new_vertices[:,0],new_vertices[:,2]-new_vertices[:,0])
        guard=survives&(np.linalg.norm(old_norm,axis=1)>1e-12)
        reversed_faces+=int(np.sum(np.einsum('ij,ij->i',old_norm[guard],new_norm[guard])<=0))
        V[a]=position;F[F==b]=a;members[a]=group
    keep=np.array([len(set(f))==3 for f in F]);F=F[keep];used=np.unique(F)
    mapping=np.full(len(V),-1,dtype=np.int64);mapping[used]=np.arange(len(used))
    matches=np.array_equal(V[used],output_V) and np.array_equal(mapping[F],output_F)
    if not matches or reversed_faces or worst>history['tolerance_mm']+1e-10:errors.append('Independent edge-collapse replay failed')
    return dict(replay_matches_exact_export=matches,maximum_original_cluster_displacement_mm=worst,
        surviving_guarded_face_reversals=reversed_faces,operations=len(history['operations']))


loader=[]
for module in ['infraspinatus','subscapularis']:
    row=parent['sources']['tpu:'+module]
    solid,V,F,meta=B.load_boundary(row)
    soup,faces=mesh_geom.read_stl(ROOT/row['file']);reference_v,reference_f=mesh_geom.weld(soup,faces,decimals=4)
    same=np.array_equal(V,reference_v) and np.array_equal(F,reference_f)
    if not same:errors.append('Legacy STL input changed')
    loader.append(dict(module=module,legacy_STL_vertices_faces_unchanged=same,source_sha256=row['sha256']))

cases=[]
for branch,subfolder in [('s1_tpu_edge_collapse_revision','epsilon_0.001'),('s1_tpu_local_fairing_revision','move_0.02')]:
    for module in ['infraspinatus','subscapularis']:
        folder=A.DESIGN/branch/module/subfolder;file=folder/'surface_revision.json';revision=read(file)
        if revision['parent_sha256']!=sha(assembly_file):errors.append('Stale surface revision')
        script=ROOT/'mujoco_arm/tools'/('build_shoulder_edge_collapse_surface.py' if 'edge_collapse' in branch else 'build_shoulder_local_fairing_surface.py')
        if sha(script)!=revision['script_sha256']:errors.append('Changed surface generator')
        V,F,actual=raw_indexed(revision['output']);old=S['tpu:'+module]
        if 'edge_collapse' in branch:
            independent=replay_collapse(old,revision['collapse_history'],V,F)
        else:
            iv,iff,_=raw_indexed(revision['input_mesh']);distance=np.linalg.norm(V-iv,axis=1)
            a0=iv[iff];a1=V[F];n0=np.cross(a0[:,1]-a0[:,0],a0[:,2]-a0[:,0]);n1=np.cross(a1[:,1]-a1[:,0],a1[:,2]-a1[:,0])
            reversed_count=int(np.sum(np.einsum('ij,ij->i',n0,n1)<0))
            if not np.array_equal(F,iff) or distance.max()>revision['maximum_allowed_vertex_displacement_mm']+1e-10:errors.append('Fairing bound/index changed')
            if reversed_count!=revision['faces_with_reversed_normal_from_input']:errors.append('Fairing reversal count changed')
            independent=dict(maximum_displacement_from_edge_collapse_input_mm=float(distance.max()),
                changed_vertex_count=int(np.sum(distance>1e-12)),reversed_faces_from_input=reversed_count)
        rigid=[];soft=[]
        for key,body in S.items():
            if key=='tpu:'+module:continue
            volume=A.I.intersection_volume(actual,body)
            if volume>.001:(soft if is_compliant(key) else rigid).append(dict(part=key,volume_mm3=volume))
        wire_hits=[]
        for route in routes:
            volume=A.I.intersection_volume(actual,A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3))
            if volume>.001:wire_hits.append(dict(head=route['head'],volume_mm3=volume))
        route=next(r for r in routes if r['module']==module);guards=[]
        for port in route['ports'].values():
            c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source']);guards.append(rod(c-n*8.,c+n*3.,12.))
        for path in [route['core_path_source_mm']]+[route['connectors'][side]['path_source_mm'] for side in ['origin','distal']]:
            guards.append(A.D.outer_envelope(np.array(path),1.8))
        guard=A.I._manifold().Manifold.batch_boolean(guards,A.I._manifold().OpType.Add)
        protected=old^guard
        added=float((actual-old).volume());removed=float((old-actual).volume())
        if max(abs(added-revision['added_material_signed_mm3']),abs(removed-revision['removed_material_signed_mm3']))>.001:errors.append('Material difference replay mismatch')
        lumen=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.86)
        cases.append(dict(branch=branch,module=module,revision_report=file.relative_to(ROOT).as_posix(),revision_report_sha256=sha(file),
            mesh_sha256=revision['output']['sha256'],independent=independent,triangles=len(F),
            added_material_signed_mm3=added,removed_original_material_signed_mm3=removed,
            protected_original_stock_loss_signed_mm3=float((protected-actual).volume()),
            nominal_lumen_overlap_before_mm3=A.I.intersection_volume(old,lumen),
            nominal_lumen_overlap_after_mm3=A.I.intersection_volume(actual,lumen),
            rigid_objects_checked=sum(not is_compliant(k) for k in S),rigid_hits_over_threshold=rigid,
            undeformed_soft_hits_over_threshold=soft,active_wire_hits_over_threshold=wire_hits,
            positive_components=sum(p.volume()>1e-6 for p in actual.decompose()),adopted_as_current=False))
        print('Verified surface',branch,module,'rigid',len(rigid),'wires',len(wire_hits),flush=True)

capture=read(HERE/'captured_attempts.json');attempts=[]
if not capture['complete'] or len(capture['runs'])!=6:errors.append('Incomplete native captures')
for row in capture['runs']:
    attempt_file=ROOT/row['attempt'];data=read(attempt_file)
    if sha(attempt_file)!=row['attempt_sha256'] or sha(ROOT/row['engine_log'])!=row['engine_log_sha256']:errors.append('Captured attempt/log changed')
    snapshot=attempt_file.parent/('generator_snapshot_'+data['script_sha256'][:12]+'.py')
    if sha(snapshot)!=data['script_sha256']:errors.append('Mesher snapshot missing')
    if row['process_exit_code']!=2 or data['geometry_tetrahedralization_pass'] or data['numerical_elasticity_mesh_quality_pass']:errors.append('Unexpected mesh state requires review')
    V,F,_=raw_indexed(data['source']);loaded=np.load(attempt_file.parent/'input_boundary.npz')
    if not np.array_equal(V,loaded['nodes_local_mm']) or not np.array_equal(F,loaded['triangles']):errors.append('Mesher changed indexed boundary')
    if data['input_preprocessing']['validation_weld_key_decimals'] is not None or data['input_preprocessing']['dropped_degenerate_triangles']!=0:errors.append('Unexpected indexed preprocessing')
    attempts.append(dict(**row,error=data['error'],error_type=data['error_type'],input_indexed_boundary_preserved=True))

remesh_file=A.DESIGN/'s1_tpu_local_fairing_revision/infraspinatus/move_0.02/surface_remesh/infraspinatus/angle_60_h0.6/mesh_attempt.json'
remesh=read(remesh_file)
if remesh['geometry_tetrahedralization_pass'] or remesh['numerical_elasticity_mesh_quality_pass']:errors.append('Unexpected remesh success')
if sha(ROOT/remesh['reference_tutorial'])!=remesh['reference_tutorial_sha256']:errors.append('Reference tutorial changed')
result=dict(scope=__doc__,script_sha256=sha(Path(__file__)),selected_assembly_sha256=sha(assembly_file),selected_source_count=len(S),
    candidate_pointer_unchanged=sha(selection_file)==before['design/S1_CURRENT_CANDIDATE.json'],legacy_loader_checks=loader,
    independent_candidate_checks=cases,captured_mesh_attempts=attempts,
    remeshing=dict(file=remesh_file.relative_to(ROOT).as_posix(),sha256=sha(remesh_file),error=remesh['error'],classified_entity_counts=remesh['classified_entity_counts']),
    complete_accepted_TPU_volume_mesh=False,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,errors=errors)
(HERE/'verification_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print('Verified',len(cases),'complete surfaces,',len(attempts),'native failure captures; errors',errors,flush=True)
if errors:raise SystemExit(1)
