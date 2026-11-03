"""Independently audit the bore redesign and reject invalid tetrahedra."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'.deps'));sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
import build_shoulder_tpu as T
from build_shoulder_shaft_collar import rod
from check_shoulder_front_liner_keepers import is_compliant


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

read=lambda p:json.loads(p.read_text(encoding='utf8'))
folder=A.DESIGN/'s1_tpu_channel_rebuild_revision/infraspinatus/radius_0.9'
file=folder/'channel_revision.json';revision=read(file);parent=read(ROOT/revision['parent_report']);errors=[]
source=revision['source'];output=revision['output']
if sha(ROOT/output['file'])!=output['sha256'] or sha(ROOT/source['file'])!=source['sha256']:raise ValueError('Source changed')
old=A.load_checked(ROOT/source['file'],source['translation_source_mm'],source['sha256'])[0]
V,F=mesh_geom.load_obj(str(ROOT/output['file']));_,boundary,nm=mesh_geom.edge_census(F)
actual=A.I._manifold().Manifold(A.I._manifold().Mesh64(np.ascontiguousarray(V),np.ascontiguousarray(F,dtype=np.uint64)))
routes=read(ROOT/revision['routes'])['routes'];route=next(r for r in routes if r['module']=='infraspinatus')
if sha(ROOT/revision['routes'])!=revision['routes_sha256']:raise ValueError('Changed routes')
path=np.array(route['full_wire_path_source_mm']);radius=revision['new_cylinder_diameter_mm']/2;node_radius=revision['new_joint_sphere_diameter_mm']/2
parts=[T.rod(a,b,radius,segments=32) for a,b in zip(path,path[1:]) if np.linalg.norm(b-a)>1e-6]
parts +=[A.I._manifold().Manifold.sphere(node_radius,32).translate(p) for p in path]
cutter=A.I._manifold().Manifold.batch_boolean(parts,A.I._manifold().OpType.Add)
expected=old-cutter;extra=float((actual-expected).volume());lost=float((expected-actual).volume())
if max(extra,lost)>.001 or boundary or nm:errors.append('Independent bore/topology replay failed')
ends=[]
for port in route['ports'].values():
    c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source']);ends.append(rod(c-n*8,c+n*3,12.))
end_guard=A.I._manifold().Manifold.batch_boolean(ends,A.I._manifold().OpType.Add)
paths=[route['core_path_source_mm']]+[route['connectors'][side]['path_source_mm'] for side in ['origin','distal']]
wall=A.I._manifold().Manifold.batch_boolean([A.D.outer_envelope(np.array(p),node_radius+.75) for p in paths],A.I._manifold().OpType.Add)-cutter-end_guard
missing=float((wall-actual).volume());outside=float(((old-actual)-cutter).volume())
if missing>.001 or outside>.001:errors.append('Bore wall/protected removal replay failed')
# Circumscribe the guard primitives so polygon inscribing cannot silently
# underestimate the requested Euclidean radial envelope.
unit=A.I._manifold().Manifold.sphere(1.,32)
sv,sf=mesh_geom.manifold_mesh_to_arrays(unit.to_mesh64());tri=sv[sf]
normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
minimum_plane_distance=float(np.min(np.abs(np.einsum('ij,ij->i',normal,tri[:,0]))))
guard_radius=node_radius+.75;guard_parts=[]
for p in paths:
    p=np.array(p)
    guard_parts +=[T.rod(a,b,guard_radius/np.cos(np.pi/32),segments=32) for a,b in zip(p,p[1:]) if np.linalg.norm(b-a)>1e-6]
    guard_parts +=[unit.scale([guard_radius/minimum_plane_distance]*3).translate(point) for point in p]
conservative=A.I._manifold().Manifold.batch_boolean(guard_parts,A.I._manifold().OpType.Add)-cutter-end_guard
conservative_missing=float((conservative-actual).volume())
if conservative_missing>.001:errors.append('Circumscribed wall guard is missing material')
rigid=[];soft=[];count=0
for key,row in parent['sources'].items():
    if key=='tpu:infraspinatus':continue
    part=A.load_checked(ROOT/row['file'],row['translation_source_mm'],row['sha256'])[0]
    volume=A.I.intersection_volume(part,actual)
    if not is_compliant(key):count+=1
    if volume>.001:(soft if is_compliant(key) else rigid).append(dict(part=key,volume_mm3=volume))
wires=[]
for route_row in routes:
    volume=A.I.intersection_volume(actual,A.D.outer_envelope(np.array(route_row['full_wire_path_source_mm']),.3))
    if volume>.001:wires.append(dict(head=route_row['head'],volume_mm3=volume))
if rigid or wires:errors.append('Rebuilt bore interferes with rigid parts/wires')

dest=folder/'volume_mesh/infraspinatus/preserve_h0.8_facetangle1e-06';attempt_file=dest/'mesh_attempt.json';attempt=read(attempt_file)
raw=np.load(dest/'raw_engine_connectivity.npz');mesh=np.load(dest/'volume_mesh.npz')
if sha(dest/'raw_engine_connectivity.npz')!=attempt['raw_engine_connectivity_sha256']:errors.append('Changed raw engine data')
lookup={int(tag):point for tag,point in zip(raw['node_tags'],raw['nodes_local_mm'])}
points=np.array([[lookup[int(tag)] for tag in cell] for cell in raw['tetrahedron_node_tags']])
matrix=np.stack([points[:,1]-points[:,0],points[:,2]-points[:,0],points[:,3]-points[:,0]],axis=2)
volumes=np.linalg.det(matrix)/6.;stored=mesh['signed_volumes_mm3']
if not np.allclose(volumes,stored,rtol=1e-8,atol=1e-8):errors.append('Independent tag lookup/signed volumes differ')
edge2=sum(np.sum((points[:,a]-points[:,b])**2,axis=1) for a in range(4) for b in range(a))
quality=12*(3*np.abs(volumes))**(2/3)/edge2
bad=int(np.sum(volumes<=0));low=int(np.sum(quality<.05))
if bad!=attempt['metrics']['negative_or_zero_volume_count'] or low!=attempt['metrics']['mean_ratio_below_0_05_count']:errors.append('Independent mesh rejection differs')
if not bad or not low or attempt['geometry_tetrahedralization_pass'] or attempt['numerical_elasticity_mesh_quality_pass']:errors.append('Invalid mesh promoted to pass')
surface=V[F];surface_volume=float(np.einsum('ij,ij->i',surface[:,0],np.cross(surface[:,1],surface[:,2])).sum()/6.)
tetra=mesh['tetrahedra'];facets=np.vstack([tetra[:,[0,1,2]],tetra[:,[0,1,3]],tetra[:,[0,2,3]],tetra[:,[1,2,3]]])
_,counts=np.unique(np.sort(facets,axis=1),axis=0,return_counts=True)
diagnostic=read(folder/'readback_surface_diagnostic.json')
if diagnostic['source_revision_sha256']!=sha(file):errors.append('Stale surface diagnostic')
if diagnostic['surface_audit']['count']!=0:errors.append('Unexpected new diagnostic crossings')

crash_dir=folder/'volume_mesh/infraspinatus/preserve_h0.8_facetangle1e-06_algorithm10'
crash=read(crash_dir/'process_result.json')
if sha(crash_dir/'engine_stdout.txt')!=crash['engine_log_sha256'] or crash['accepted_volume_mesh']:errors.append('Crash evidence changed')
abort_dir=folder/'surface_remesh/infraspinatus/angle_60_h0.6';abort=read(abort_dir/'process_result.json')
if sha(abort_dir/'engine_stdout.txt')!=abort['engine_log_sha256'] or abort['accepted_mesh']:errors.append('Interrupted remesh evidence changed')
guard=read(ROOT/'reviews/S1_2026-11-03/indexed_boundary_meshing/budget_control/supervised/process_result.json')
if guard['stop_reason']!='native_stdout_resource_limit' or not guard['worker_terminal']:errors.append('Resource control did not reap worker')
result=dict(scope=__doc__,script_sha256=sha(Path(__file__)),source_revision=file.relative_to(ROOT).as_posix(),source_revision_sha256=sha(file),
    indexed_export_boundary_edges=boundary,indexed_export_nonmanifold_edges=nm,
    independent_bore_replay_extra_mm3=extra,independent_bore_replay_loss_mm3=lost,
    explicit_removed_material_mm3=float((old-actual).volume()),unintended_loss_outside_new_cutter_mm3=outside,
    missing_nominal_0_75mm_wall_guard_outside_end_interfaces_mm3=missing,
    circumscribed_0_75mm_wall_guard=dict(missing_material_outside_end_interfaces_mm3=conservative_missing,
        cylinder_segments=32,cylinder_circumradius_mm=float(guard_radius/np.cos(np.pi/32)),
        sphere_unit_minimum_face_plane_distance=minimum_plane_distance,
        interpretation='Conservative polygonal guard around internal path; nominal envelope audit, not printing or deformed-wall proof.'),
    rigid_objects_checked=count,rigid_hits_over_threshold=rigid,undeformed_soft_hits_over_threshold=soft,active_wire_hits_over_threshold=wires,
    surface_diagnostic_detected_crossings=0,surface_diagnostic_is_not_exact_arithmetic_proof=True,
    source_coordinate_boundary_signed_volume_mm3=surface_volume,
    independent_tetra_audit=dict(attempt_file=attempt_file.relative_to(ROOT).as_posix(),attempt_sha256=sha(attempt_file),
        raw_tag_dictionary_lookup_verified=True,tetra_count=len(volumes),negative_or_zero_volumes=bad,
        minimum_signed_volume_mm3=float(volumes.min()),total_signed_volume_mm3=float(volumes.sum()),
        minimum_mean_ratio=float(quality.min()),mean_ratio_below_0_05_count=low,
        tetra_boundary_faces=int(np.sum(counts==1)),tetra_faces_with_more_than_two_uses=int(np.sum(counts>2)),
        source_volume_difference_mm3=float(volumes.sum()-actual.volume()),accepted=False),
    native_crash=crash,interrupted_reparametrization=abort,resource_budget_control_passed=True,
    adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,errors=errors)
(folder/'independent_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print('Bore replay verified; rigid/wire hits',len(rigid),len(wires),'wall missing',missing,flush=True)
print('Complete tetra mesh rejected:',bad,'nonpositive elements,',low,'below quality 0.05; errors',errors,flush=True)
if errors:raise SystemExit(1)
