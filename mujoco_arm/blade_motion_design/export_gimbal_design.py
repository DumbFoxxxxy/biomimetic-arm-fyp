"""Export the requested no-payload motion baseline with explicit hardware limits."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import copy,csv,json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from gimbal_arm_model import build_gimbal,GH_SHIFT
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from verification_snapshot import fingerprint
from shoulder_model import EXPORT,v8

def main():
    tree,specs=build_gimbal();root=tree.getroot();current=fingerprint(root);reports={}
    for name in ['gimbal_assembly_report','gimbal_wire_geometry_report','gimbal_wire_pair_report','gimbal_inverse_drive_report']:
        report=json.loads((EXPORT/f'{name}.json').read_text())
        if not report['pass'] or report.get('input_fingerprint')!=current:raise RuntimeError(f'{name}: failed or stale; do not freeze anchors')
        reports[name]=report
    gravity_tree=copy.deepcopy(root);gravity_tree.find('option').set('gravity','0 0 -9.81');gravity=json.loads((EXPORT/'gimbal_gravity_drive_report.json').read_text())
    if gravity.get('input_fingerprint')!=fingerprint(gravity_tree):raise RuntimeError('Gravity diagnostic is stale')
    for node in gravity_tree.find('actuator'):node.set('forcerange','-300 0')
    conditional_path=EXPORT/'gimbal_gravity_300N_drive_report.json';conditional=json.loads(conditional_path.read_text()) if conditional_path.exists() else {}
    conditional_pass=bool(conditional.get('pass') and conditional.get('input_fingerprint')==fingerprint(gravity_tree))
    connectivity=json.loads((EXPORT/'gimbal_connectivity_report.json').read_text())
    if any(r['positive_material_component_count']!=1 for r in connectivity):raise RuntimeError('Disconnected carrier material')
    m=_load_project_xml(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);targets=target_actions();zero=next(iter(targets.values()))[0];apply_target(m,d,zero)
    ghid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'gh_flex');gh_world=v8.mj_to_design(d.xpos[ghid]);names=[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,t) for t in range(m.ntendon)]
    rows=[]
    for tid,rid in enumerate(names):
        node=root.find(f'./tendon/spatial[@name="{rid}"]');first=node.find('site');last=list(node)[-1];sn=first.get('site');en=last.get('site');a=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,sn);b=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,en)
        aa=v8.mj_to_design(d.site_xpos[a]);bb=v8.mj_to_design(d.site_xpos[b]);ab=int(m.site_bodyid[a]);eb=int(m.site_bodyid[b])
        al=v8.mj_to_design(m.site_pos[a]);bl=v8.mj_to_design(m.site_pos[b])
        row={'id':rid,'supply_exit_site':sn,'supply_carrier':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,ab) or 'world','termination_site':en,'termination_carrier':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,eb),'neutral_active_length_mm':float(d.ten_length[tid]*1000)}
        for label,p in [('exit_zero_world',aa),('end_zero_world',bb),('exit_body_local',al),('end_body_local',bl)]:
            for axis,value in zip('xyz',p):row[f'{label}_{axis}_mm']=float(value)
        rows.append(row)
    lengths=[];render=[];rom=[];hid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'humerus');uid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'ulna');apply_target(m,d,zero);neutral=(d.xpos[uid]-d.xpos[hid]);neutral/=np.linalg.norm(neutral)
    for action,seq in targets.items():
        dense=[a*(1-f)+b*f for a,b in zip(seq[:-1],seq[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[seq[-1]]
        for z in dense:apply_target(m,d,z);lengths.append(d.ten_length.copy()*1000)
        direction=d.xpos[uid]-d.xpos[hid];direction/=np.linalg.norm(direction);rom.append({'action':action,'endpoint_parameters_rad':seq[-1].tolist(),'world_arm_angle_from_zero_deg':float(np.degrees(np.arccos(np.clip(direction@neutral,-1,1))))})
        if action in ['SCAP_UP_ROT','GH_FLEX','GH_ABD','ELBOW','RAISE_170']:render.append({'action':action,'qpos':d.qpos.tolist()})
    lengths=np.array(lengths)
    for i,row in enumerate(rows):row.update(minimum_active_length_mm=float(lengths[:,i].min()),maximum_active_length_mm=float(lengths[:,i].max()),required_geometric_excursion_mm=float(np.ptp(lengths[:,i])))
    with (EXPORT/'full_arm_cable_coordinates.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (EXPORT/'gimbal_render_states.json').write_text(json.dumps(render,indent=2));(EXPORT/'full_arm_rom.json').write_text(json.dumps(rom,indent=2));tree.write(EXPORT/'gimbal_arm.xml',encoding='utf-8',xml_declaration=True)
    bundle={'scope':'Verified sampled motion-design baseline under the user no-load acceptance scope; not a manufacturing release or a built prototype','motion_design_gate_pass':True,'hardware_ready':False,'coordinate_units':'mm; x lateral, y superior, z posterior; body-local axes use the same coordinate convention before body rotation','gh_shift_design_mm':GH_SHIFT.tolist(),'gh_zero_world_mm':gh_world.tolist(),'scalar_dofs':m.nv,'tendon_count':m.ntendon,'shoulder_tendon_count':sum(n.startswith('S') for n in names),'steel_diameter_mm':1.,'terminal_contact_zone_mm':5.,'axes':specs,'motion_count':len(targets),'geometry_pose_samples':reports['gimbal_assembly_report']['exact_pose_checkpoints'],'wire_pair_gap_mm':reports['gimbal_wire_pair_report']['minimum_polyline_surface_pair_gap_mm'],'estimated_self_weight_200N_diagnostic_pass':gravity['pass'],'estimated_self_weight_failed_actions':[r['action'] for r in gravity['actions'] if not r['pass']],'gravity_diagnostic':'Uniform solid printed material density estimate; no payload; recorded separately from requested motion gates; not actuator sizing or strength qualification','input_fingerprint':current}
    (EXPORT/'full_arm_motion_design_checkpoint.json').write_text(json.dumps(bundle,indent=2));print({k:v for k,v in bundle.items() if k not in ['axes','input_fingerprint']})
    bundle['estimated_self_weight_300N_capability_diagnostic_pass']=conditional_pass
    bundle.update(actuation_scope='Ideal pull-only tendons; TPU muscle structure and deformation are not modeled',user_intended_tpu_muscle_actuation_verified=False,full_fyp_solution=False,scope_correction='design/TPU_MUSCLE_SCOPE_CORRECTION.md')
    if conditional_pass:
        bundle['estimated_self_weight_peak_tension_N']=max(r['maximum_feedforward_tension_N'] for r in conditional['actions'])
        ET.ElementTree(gravity_tree).write(EXPORT/'gimbal_gravity_300N.xml',encoding='utf-8',xml_declaration=True)
    (EXPORT/'full_arm_motion_design_checkpoint.json').write_text(json.dumps(bundle,indent=2))
if __name__=='__main__':main()
