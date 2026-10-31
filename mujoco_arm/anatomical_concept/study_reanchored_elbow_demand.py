"""Actual endpoints, chord demand and chord/bone diagnostics, not muscle strain.

No deformation law or virtual guide is imposed. These diagnostics decide which
poses still need real compliant contact/wrapping analysis before drive claims.
"""
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation

def run():
    assert EXTENSOR_REANCHOR
    a=json.loads((DEST/'assembly.json').read_text());j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);rows={r['key']:r for r in a['instances']}
    # Only actual bone solids: a chord is not the actual wire inside a curved
    # belly and its end frames. Do not collide it with its own end hardware.
    bones={k:load(ROOT/rows[k]['file'])[0].translate(rows[k]['translation_display_mm']) for k in ['FJ3368','FJ3391']}
    modules=[]
    for b in a['muscle_end_boundaries']:
        A=np.array(b['proximal_pad_centre_display_mm']);B=np.array(b['distal_pad_zero_centre_display_mm']);samples=[]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];D=E+R@(B-E);L=float(np.linalg.norm(A-D));moment=float(np.cross(D-E,unit(A-D))@axis);wire=rod(A,D,.3);hits=[]
            for key,s in bones.items():
                target=s.transform(T) if rows[key]['owner']=='elbow' else s
                overlap=intersection(wire,target)
                if overlap>.001:hits.append({'bone':key,'volume_mm3':overlap})
            samples.append({'elbow_deg':angle,'distal_pad_display_mm':D.tolist(),'chord_mm':L,'chord_moment_mm':moment,'straight_chord_bone_hits':hits})
        lengths=[s['chord_mm'] for s in samples];modules.append({'id':b['muscle'],'zero_chord_mm':lengths[0],'min_chord_mm':min(lengths),'max_chord_mm':max(lengths),'max_endpoint_chord_shortening_fraction':1-min(lengths)/lengths[0],'max_endpoint_chord_growth_fraction':max(lengths)/lengths[0]-1,'minimum_correct_chord_moment_mm':min((1 if b['muscle']=='FJ1486' else -1)*s['chord_moment_mm'] for s in samples),'first_straight_chord_bone_hit_deg':next((s['elbow_deg'] for s in samples if s['straight_chord_bone_hits']),None),'poses':samples})
        print(b['muscle'],'chordchange',max(lengths)-min(lengths),'firstbonehit',modules[-1]['first_straight_chord_bone_hit_deg'],flush=True)
    (DEST/'endpoint_demand_diagnostics.json').write_text(json.dumps({'scope':'endpoint chord and straight-line bone diagnostics only; not actual steel travel, TPU strain, compliant contact or drive validation','angle_step_deg':1,'modules':modules,'frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
