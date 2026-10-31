"""Source-guarded pin retention candidate, explicitly not an actuated elbow pass."""
from build_direct_muscle_yoke import *
from pathlib import Path

def run():
    dest=DEST/'paired_swivel_support/captive_keepers';read=lambda n:json.loads((dest/n).read_text(encoding='utf8'));digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();a=read('assembly.json');d=read('design.json');g=read('geometry_checks.json');fit=read('assembly_paths.json');c=read('cad_roundtrip.json');whole=read('whole_arm_STEP_check.json');mesh=read('export_checks.json');assert len(a['instances'])==127 and d['selected'] is not None
    assert not g['neutral_rigid_interferences'] and not g['neutral_TPU_contact_flags'] and g['first_return_interference'] is None and g['motion_attached_rigid_count']==123 and g['return_samples']==266
    for r in [d,g,fit]:
        for path,h in r['source_sha256'].items():assert digest(Path(path))==h,(path,'stale source')
    assert len(d['capture_witnesses'])==4
    for witness in d['capture_witnesses']:assert witness['rest_keeper_overlap_mm3']<.001 and witness['outward_0_09mm_keeper_overlap_mm3']<.001 and witness['outward_0_11mm_keeper_overlap_mm3']>.001 and witness['inward_0_01mm_washer_overlap_mm3']>.001 and witness['minimum_stub_engagement_at_outward_stop_mm']>=3.2
    assert len(fit['paths'])==19 and sum(p['rigid_clear'] for p in fit['paths'])==15 and fit['bare_forearm_mount_found']
    assert len(c['parts'])==len(mesh['parts'])==15 and all(p['pass'] for p in c['parts']+mesh['parts'])
    assert whole['valid'] and whole['solid_count']==125 and whole['source_assembly_sha256']==mesh['assembly_sha256']==digest(dest/'assembly.json') and whole['STEP_sha256']==digest(ROOT/whole['STEP'])
    for r in whole['parts'].values():assert r['OBJ_sha256']==digest(ROOT/r['source_OBJ']) and r['STEP_sha256']==digest(ROOT/r['source_STEP']) and r['STEP_vs_OBJ_vertex_error_mm']<.001
    for scene,count in [('neutral_view.xml',127),('rigid_hinges.xml',123)]:
        for pose in mesh['compiled_scenes'][scene]:assert pose['registered_count']==count and pose['maximum_error_mm']<.001
    names=['assembly.json','design.json','geometry_checks.json','assembly_paths.json','cad_roundtrip.json','whole_arm_STEP_check.json','export_checks.json','rigid_model_scope.json','render_inputs.json'];result={'stage':'M10-R1 Ukeeper captive pivot candidate','nominal_neutral_meshes':127,'known_rigid_motion_parts':123,'nominal_neutral_STEP_solids':125,'new_revised_parts':15,'conditional_forward_and_return_samples':266,'selected_rigid_subassembly_paths':15,'rejected_access_trials':4,'pin_capture_nominal_outward_float_mm':.1,'minimum_stub_engagement_at_capture_mm':3.2,'positive_capture_geometry_checked':True,'physical_retention_or_loosening_tested':False,'thin_sheet_fabrication_verified':False,'complete_arm_or_TPU_assembly_verified':False,'TPU_deformation_verified':False,'actual_internal_steel_motion_verified':False,'actual_actuated_elbow_range_deg':None,'frozen':False,'evidence_sha256':{n:digest(dest/n) for n in names}};(dest/'review_status.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('captive keeper geometry/CAD/subassembly review passed, no force/drive claim',flush=True)

if __name__=='__main__':run()
