"""Compare a direct bone-to-TPU-end interface before adding more cable hardware.

Only endpoint kinematics and a straight chord force approximation. No lattice
deformation, bone clearance, actual wire length or actuation acceptance.
"""
from build_forearm_receivers import *

def run():
    folder=OUT/'direct_muscle_interfaces';folder.mkdir(exist_ok=True);a=json.loads((PARENT/'assembly.json').read_text());offsets={p['id']:np.array(p['translation_display_mm']) for p in a['placements']};j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);reports=[]
    for fid in ['FJ1486','FJ1477']:
        source=OUT/'elbow_muscles'/fid/'design.json';d=json.loads(source.read_text());A,B=np.array(d['pad_centres_display_mm'])+offsets[fid];length0=np.linalg.norm(B-A);rows=[]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();D=E+R@(B-E);length=np.linalg.norm(A-D);force=unit(A-D);moment=float(np.cross(D-E,force)@axis);rows.append({'elbow_angle_deg':angle,'distal_pad_centre_display_mm':D.tolist(),'endpoint_chord_mm':float(length),'chord_change_from_zero_mm':float(length-length0),'chord_force_moment_arm_mm':moment})
        correct=(lambda m:m>0) if fid=='FJ1486' else (lambda m:m<0);firstwrong=next((r['elbow_angle_deg'] for r in rows if not correct(r['chord_force_moment_arm_mm'])),None);reports.append({'muscle':fid,'proximal_pad_display_mm':A.tolist(),'distal_pad_zero_display_mm':B.tolist(),'zero_chord_mm':float(length0),'minimum_chord_mm':min(r['endpoint_chord_mm'] for r in rows),'maximum_chord_mm':max(r['endpoint_chord_mm'] for r in rows),'largest_endpoint_chord_shortening_fraction':float((length0-min(r['endpoint_chord_mm'] for r in rows))/length0),'largest_endpoint_chord_lengthening_fraction':float((max(r['endpoint_chord_mm'] for r in rows)-length0)/length0),'first_wrong_chord_force_sign_deg':firstwrong,'poses':rows,'source_design_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
    result={'scope':'proposed direct mechanical connection of distal TPU end to forearm, pad position rotates with forearm; no such adapter constructed yet. Chord metrics are not muscle strain or steel excursion. Straight chord force sign cannot replace actual curved TPU/contact analysis.','modules':reports,'bone_hard_parts_changed':False,'endpoint_constraint_applied_to_current_model':False,'muscle_deformation_solved':False,'drive_verified':False,'frozen':False};(folder/'endpoint_chord_study.json').write_text(json.dumps(result,indent=2),encoding='utf8');print([{k:r[k] for k in ['muscle','zero_chord_mm','largest_endpoint_chord_shortening_fraction','largest_endpoint_chord_lengthening_fraction','first_wrong_chord_force_sign_deg']} for r in reports],flush=True)

if __name__=='__main__':run()
