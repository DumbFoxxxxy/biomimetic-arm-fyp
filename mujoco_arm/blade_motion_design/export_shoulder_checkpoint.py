"""Export the verified shoulder subset without freezing full-arm hardpoints."""
import csv,json
from shoulder_model import inputs,EXPORT,SOURCE,input_fingerprints

def export():
    cfg,fit,anchors=inputs();geometry=json.loads((EXPORT/'shoulder_geometry_report.json').read_text());drive=json.loads((EXPORT/'shoulder_drive_report.json').read_text());channel=json.loads((EXPORT/'channel_clearance_report.json').read_text())
    expected=input_fingerprints()
    for report in [geometry,drive,channel]:
        if report['input_sha256']!=expected:raise RuntimeError('Report/input hash mismatch')
    okay=geometry['pass'] and drive['pass'] and channel['free_channel_pass'] and channel['cable_cable_pass']
    with (EXPORT/'shoulder_anchor_checkpoint.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['scope','cable','part','point','frame','x_mm','y_mm','z_mm','minimum_length_mm','maximum_length_mm','excursion_mm'])
        for r in anchors['routes']:
            for part,kind,frame,key in [('thorax_frame','fixed_eye','world_design','fixed_eye_world_mm'),('scapula','moving_eye','fitted_zero_world','moving_eye_zero_world_mm'),('scapula','moving_eye','AC_local_design','eye_local_ac_mm')]:
                writer.writerow(['SHOULDER_SUBSYSTEM_ONLY',r['id'],part,kind,frame,*r[key],r['min_length_mm'],r['max_length_mm'],r['required_excursion_mm']])
    result={'shoulder_subsystem_pass':bool(okay),'full_arm_design_complete':False,'frozen_full_arm_hardpoints':False,
        'reference':'A Shoulder Structure of Muscle-Driven Humanoid with Shoulder Blades, local PDF pp.3-5',
        'geometry_change':{'scapula_outward_mm':12.,'clavicle_anterior_mm':8.,'joint_centers_preserved':True,'local_channel':'S07'},
        'nominal_cable_diameter_mm':1.,'terminal_eye_clamp_contact_length_mm':3.,'channel_sweep_radius_mm':2.,
        'shoulder_motion_targets':{'upward_rotation_deg':50.,'downward_rotation_deg':20.,'protraction_retraction_mm':15.,'elevation_depression_mm':20.},
        'maximum_tracking_error_mm':max(r['maximum_contact_tracking_error_mm'] for r in drive['actions']),
        'maximum_tracking_error_deg':max(r['maximum_orientation_tracking_error_deg'] for r in drive['actions']),
        'channel_wall_clearance_mm':channel['minimum_channel_wall_clearance_mm'],'cable_surface_pair_gap_mm':channel['minimum_cable_cable_surface_gap_mm'],
        'pending':['GH cable anchors and tendon-driven ROM','limb/scapula and limb/cable interactions','coupled adduction and combined reaching','final design and manufacturing export'],
        'input_sha256':expected}
    (EXPORT/'shoulder_checkpoint.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':export()
