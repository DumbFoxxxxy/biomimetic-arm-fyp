"""Elbow handoff evidence; deliberately does not certify compliant actuation."""
from build_direct_muscle_yoke import *

def read(path):return json.loads(path.read_text(encoding='utf8'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    assert ('--access-relief' in sys.argv or EXTENSOR_REANCHOR) and SHIFTED
    a=read(DEST/'assembly.json');cad=read(DEST/'whole_arm_STEP_check.json');audit=read(DEST/'bone_motion_audit.json');paths=read(DEST/'motion_paths.json');neutral=read(DEST/'interface_checks.json');mesh=read(DEST/'mesh_export_checks.json');fit=read(DEST/'assembly_checks.json')
    assert len(a['instances'])==107 and cad['valid'] and cad['solid_count']==105
    assert cad['source_assembly_sha256']==digest(DEST/'assembly.json')
    assert cad['STEP_sha256']==digest(ROOT/cad['STEP'])
    for r in cad['parts'].values():
        assert digest(ROOT/r['source_OBJ'])==r['OBJ_sha256']
        assert digest(ROOT/r['source_STEP'])==r['STEP_sha256']
        assert r['STEP_vs_OBJ_vertex_error_mm']<.001
    assert audit['complete_prescribed_suite'] and len(paths)==13
    for name,h in audit['input_evidence_sha256'].items():assert digest(DEST/name)==h
    maxstep=0
    for key,path in paths.items():
        q=np.array([r['qpos'] for r in path]);action=next(r for r in audit['actions'] if r['probe']==key)
        assert action['forward_samples']==action['return_samples']==len(path)
        for start in [0,4,8]:
            sub=q[:,start:start+4];assert abs(np.linalg.norm(sub,axis=1)-1).max()<1e-8
            maxstep=max(maxstep,float(np.degrees(2*np.arccos(np.clip(abs(np.sum(sub[:-1]*sub[1:],axis=1)),0,1))).max()))
        maxstep=max(maxstep,float(abs(np.degrees(np.diff(q[:,-1]))).max()))
        requested=np.array(action['requested_controls_deg']);index=np.flatnonzero(abs(requested)>1e-8)[0]
        for pose in path:
            if pose['controls_deg'][index]/requested[index]<=action['continuous_prefix_fraction']+1e-8:assert not pose['hits']
    assert maxstep<=1.000001
    elbow=next(r for r in audit['actions'] if r['probe']=='elbow_flexion');assert elbow['continuous_prefix_fraction']==1
    assert not neutral['neutral_rigid_interferences'] and not neutral['neutral_TPU_review_flags']
    assert all(r['nominal_contact_pass'] for r in neutral['bone_mount_seats'])
    assert all(not r['nominal_threading_path_hits'] for r in neutral['steel_neutral_threading'])
    assert mesh['all_roundtrips_pass'] and not mesh['neutral_reference_intersections']
    assert len(mesh['compiled_scenes']['arm.xml'])==103 and len(mesh['compiled_scenes']['neutral_view.xml'])==107
    assert all(r['pass'] for scene in mesh['compiled_scenes'].values() for r in scene)
    assert len(fit['paths'])==fit['rigid_clear_count']==19
    for path,h in fit['source_sha256'].items():assert digest(__import__('pathlib').Path(path))==h
    muscle_records=[]
    for fid in ['FJ1486','FJ1477']:
        folder=muscle_folder(fid)
        design=read(folder/'design.json');geometry=read(folder/'geometry_checks.json');prints=read(folder/'print_candidates/print_transforms.json')
        assert geometry['all_roundtrips_pass'] and geometry['single_connected_material_each'] and geometry['nominal_steel_threading_geometry_pass']
        for r in prints['parts']:
            assert r['print_STL_roundtrip_pass'] and digest(folder/(r['part']+'.obj'))==r['source_OBJ_sha256']
            assert digest(folder/'print_candidates'/(r['part']+'.stl'))==r['print_STL_sha256']
        boundary=next(r for r in a['muscle_end_boundaries'] if r['muscle']==fid);offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));F=np.array(design['basis_columns_display']);prox=np.array(boundary['proximal_pad_centre_display_mm']);dist=np.array(boundary['distal_pad_zero_centre_display_mm']);steel=np.array(design['wire_polyline_display_mm'])+offset;steel[-1]=dist+F[:,2]*4.6
        muscle_records.append({'id':fid,'folder':str(folder.relative_to(ROOT)),'body_span_mm':design['body_span_mm'],'TPU_volume_mm3':design['TPU_material_volume_mm3'],'basis_columns':F.tolist(),'proximal_pad_display_mm':prox.tolist(),'distal_pad_zero_display_mm':dist.tolist(),'conduit_nominal_reaction_shoulder_display_mm':(prox-F[:,2]*5.15).tolist(),'conduit_counterbore_diameter_mm':2.2,'conduit_outer_diameter_starting_candidate_mm':2.0,'conduit_grip_verified':False,'steel_diameter_candidate_mm':.6,'internal_threading_polyline_display_mm':steel.tolist(),'distal_crimp_midplane_display_mm':(dist+F[:,2]*4.6).tolist(),'crimp_barrel_blank_OD_ID_length_mm':[2.5,.9,3.0],'crimp_grip_verified':False,'installed_force_path':'proximal cap/adapter to humerus; steel through TPU guides to distal crimp seat; distal integral cap/yoke to ulna; motor/conduit installation and TPU compliance require bench verification'})
    evidence=['assembly.json','whole_arm_STEP_check.json','bone_motion_audit.json','motion_paths.json','interface_checks.json','mesh_export_checks.json','assembly_checks.json','cad_roundtrip.json','design.json' if EXTENSOR_REANCHOR else 'study.json','render_inputs.json']
    report={'stage':'M9 elbow direct TPU candidate, access relief','user_stop_boundary':'elbow handoff only; no further shoulder or whole-arm muscle work','neutral_instances':107,'attached_rigid_motion_instances':103,'neutral_rigid_STEP_solids':105,'sampled_rigid_elbow_deg':[0,150],'maximum_joint_sample_step_deg':maxstep,'boundary_refinement_deg':.1,'rigid_insertion_paths':19,'soft_insertion_observations':fit['soft_flag_count'],'complete_arm_assembly_verified':False,'TPU_deformation_verified':False,'TPU_actuated_elbow_range_deg':None,'extension_force_sign_across_range_verified':False,'motor_and_conduit_installation_verified':False,'physical_test_verified':False,'frozen':False,'evidence_sha256':{n:digest(DEST/n) for n in evidence},'muscles':muscle_records}
    if EXTENSOR_REANCHOR:
        report['stage']='M9 independent extensor reattachment candidate, not full compliant drive'
        report['endpoint_demand_diagnostics']=read(DEST/'endpoint_demand_diagnostics.json')
        report['evidence_sha256']['endpoint_demand_diagnostics.json']=digest(DEST/'endpoint_demand_diagnostics.json')
    (DEST/'review_status.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    (DEST/'elbow_interfaces.json').write_text(json.dumps({'units':'mm','coordinate_system':'display = BP3D source minus [-140,-70,1100] mm; never recenter parts separately','muscles':muscle_records,'frozen':False},indent=2),encoding='utf8')
    print('elbow candidate review passed; TPU drive, full assembly and physical performance remain unverified',flush=True)

if __name__=='__main__':run()
