"""Consolidate M10 geometry evidence without certifying muscle actuation."""
from build_direct_muscle_yoke import *
from pathlib import Path

def run():
    dest=DEST/'paired_swivel_support';read=lambda name:json.loads((dest/name).read_text(encoding='utf8'));digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest();assembly=read('assembly.json');design=read('design.json');geometry=read('geometry_checks.json');exports=read('export_checks.json');cad=read('cad_roundtrip.json');whole=read('whole_arm_STEP_check.json')
    assert len(assembly['instances'])==117 and design['selected']['positive_components']==1 and design['paired_cap_pack_first_interference'] is None
    assert design['sample_count']==geometry['return_samples']==266 and geometry['maximum_elbow_and_cap_step_deg']<=1.000001
    assert not geometry['neutral_rigid_interferences'] and not geometry['neutral_TPU_contact_flags'] and geometry['first_return_interference'] is None and geometry['motion_attached_rigid_count']==113
    for report in [design,geometry]:
        for path,h in report['source_sha256'].items():assert digest(Path(path))==h,(path,'stale geometry evidence')
    assert len(cad['parts'])==len(exports['parts'])==11 and all(r['pass'] for r in cad['parts']+exports['parts'])
    assert exports['assembly_sha256']==whole['source_assembly_sha256']==digest(dest/'assembly.json')
    assert whole['valid'] and whole['solid_count']==115 and whole['STEP_sha256']==digest(ROOT/whole['STEP'])
    for r in whole['parts'].values():
        assert r['OBJ_sha256']==digest(ROOT/r['source_OBJ']) and r['STEP_sha256']==digest(ROOT/r['source_STEP']) and r['STEP_vs_OBJ_vertex_error_mm']<.001
    for scene,expected in [('neutral_view.xml',117),('rigid_hinges.xml',113)]:
        for pose in exports['compiled_scenes'][scene]:assert pose['registered_count']==expected and pose['maximum_error_mm']<.001
    names=['assembly.json','design.json','geometry_checks.json','conditional_motion.json','cad_roundtrip.json','whole_arm_STEP_check.json','export_checks.json','rigid_model_scope.json','render_inputs.json','rigid_motion.gif']
    report={'stage':'M10 paired passive distal caps, actual bone-rooted frame','neutral_meshes':117,'attached_rigid_motion_parts':113,'nominal_neutral_STEP_solids':115,'new_native_CAD_parts':11,'conditional_forward_samples':266,'independent_return_samples':266,'maximum_elbow_and_cap_step_deg':geometry['maximum_elbow_and_cap_step_deg'],'conditional_rigid_elbow_range_deg':[0,150],'actual_passive_equilibrium_verified':False,'full_assembly_paths_verified':False,'pin_threading_retention_verified':False,'actual_internal_steel_motion_verified':False,'TPU_deformation_verified':False,'physical_test_verified':False,'actuated_elbow_range_deg':None,'frozen':False,'next_work_scope':'elbow only, stop after elbow design; no shoulder/scapular muscle expansion','evidence_sha256':{name:digest(dest/name) for name in names}}
    (dest/'review_status.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('M10geometry review passed; assembly,retention,steel/TPU anddrive remainunverified',flush=True)

if __name__=='__main__':run()
