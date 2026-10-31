"""Scope-limited M8 evidence audit; deliberately leaves M7 main index intact."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_forearm_receivers import *
from project_runtime import WorkspacePath as Path
MOUNT=DEST/'rephased_anchors'

def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def unchanged(mapping):
    for p,h in mapping.items():assert digest(Path(p))==h,p

def run():
    screen=read(MOUNT/'mounted_direction_screen.json');assert screen['first_rigid_limit'] is None and not screen['neutral_contacts'] and screen['samples']==273 and screen['maximum_step_including_receiver_hinges_deg']<=1.00000001;unchanged(screen['source_sha256']);assert sum('repaired_supports' in p for p in screen['source_sha256'])==31
    local=[]
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;r=read(folder/'unit_checks.json');assert r['all_local_motion_clear'] and r['all_local_fitting_clear'];assert all(s['nominal_contact'] for s in r['nominal_seat_contacts'])
        for k,h in r['source_sha256'].items():assert digest(folder/(k+'.obj'))==h
        c=read(folder/'cad_roundtrip.json');assert len(c['parts'])==13 and all(p['pass'] for p in c['parts']);local.append(fid)
    coupled=read(MOUNT/'coupled_path_checks.json');assert coupled['all_clear'] and coupled['screen_sha256']==digest(MOUNT/'mounted_direction_screen.json')
    for r in coupled['modules']:
        for k,h in r['source_sha256'].items():assert digest(DEST/r['muscle']/(k+'.obj'))==h
    mount=read(MOUNT/'mounting_checks.json');assert len(mount['paths'])==13 and mount['all_clear'] and mount['soft_flag_count']==0;unchanged(mount['source_sha256'])
    cad=read(MOUNT/'assembly_STEP_check.json');assert cad['solid_count']==27 and cad['valid'] and cad['relative_volume_error']<.0001 and cad['direction_screen_sha256']==digest(MOUNT/'mounted_direction_screen.json')
    assert digest(MOUNT/'cad/receiving_interfaces_review.step')==cad['STEP_sha256']
    for r in cad['parts']:assert digest(ROOT/r['source_STEP'])==r['STEP_sha256'] and digest(ROOT/r['source_OBJ'])==r['OBJ_sha256'] and r['source_vs_CAD_vertex_error_mm']<.001
    exports=read(MOUNT/'export_checks.json');assert exports['all_pass'] and len(exports['roundtrips'])==len(exports['registration'])==27 and exports['source_xml_sha256']==digest(MOUNT/'whole_arm_direction_study.xml')
    render=read(MOUNT/'render_inputs.json');assert render['mesh_count']==179 and render['screen_sha256']==digest(MOUNT/'mounted_direction_screen.json');unchanged(render['source_sha256'])
    direct=read(MOUNT/'direct_departure_probe.json');unchanged(direct['source_sha256']);assert all(r['all_tangent_point_iterations_converged'] for r in direct['modules'])
    result={'scope':'M8 conditional receiving-interface geometry, CAD and local fitting, not complete tendon/TPU actuation acceptance','parts':27,'conditional_elbow_samples':screen['samples'],'maximum_joint_step_deg':screen['maximum_step_including_receiver_hinges_deg'],'reference_bones_extra':31,'all_scoped_interface_checks_pass':True,'direct_route_collision_samples':{r['muscle']:r['guide_collision_pose_count'] for r in direct['modules']},'direct_route_verified':False,'full13_arm_action_suite_verified_for_M8':False,'complete_arm_assembly_verified':False,'TPU_deformation_verified':False,'full_arm_drive_verified':False,'physical_test_verified':False,'main_candidate_replaced':False,'frozen':False};(MOUNT/'study_status.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M8 scoped geometry checks pass; direct route remains rejected; main candidate unchanged',flush=True)

if __name__=='__main__':run()
