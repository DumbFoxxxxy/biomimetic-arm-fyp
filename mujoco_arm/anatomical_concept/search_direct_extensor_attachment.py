"""Alternative to adding a wrap guide: relocate only the TPU distal attachment.

Finite source-derived endpoint search, chord-force screening plus actual distal
hardware sweep. A correct chord sign does not prove compliant muscle torque.
The actual neutral profile must be regenerated and checked after selection.
"""
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation

def run():
    extended='--extended' in sys.argv
    dest=OUT/'direct_muscle_interfaces'/('extensor_reanchor_extended' if extended else 'extensor_reanchor_study');dest.mkdir(exist_ok=True)
    source=OUT/'direct_muscle_interfaces/yoke_access_relief_study'
    a=json.loads((source/'assembly.json').read_text());rows={r['key']:r for r in a['instances']}
    fixed={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if r['kind']!='TPU' and k!='direct_muscle_yoke' and not k.startswith('FJ1477_distal_') and k!='FJ1477_active_wire_crimp_barrel'}
    placements={p['id']:np.array(p['translation_display_mm']) for p in a['placements']}
    # The parent integrated yoke is not reusable after moving the attachment.
    # Retain the flexor's real cap for cross-module swept tests.
    flexcap=OUT/'direct_muscle_interfaces/shifted_muscle/FJ1486/distal_rigid_cap.obj'
    fixed['FJ1486_distal_rigid_cap']=load(flexcap)[0].translate(placements['FJ1486']);rows['FJ1486_distal_rigid_cap']={'owner':'elbow'}
    old=json.loads((M2/'assembly.json').read_text());moving={}
    for r in old['instances']:
        if r['key'].startswith('FJ1477_distal_') or r['key']=='FJ1477_active_wire_crimp_barrel':moving[r['key']]=load(ROOT/r['file'])[0].translate(r['translation_display_mm'])
    paths=[source/'assembly.json',M2/'assembly.json',BASE/'design.json',flexcap]
    paths += [ROOT/r['file'] for r in a['instances']]
    paths += [ROOT/r['file'] for r in old['instances'] if r['key'] in moving]
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in fixed:
            p=OUT/'repaired_supports'/(r['id']+'.obj');key='reference_'+r['id'];fixed[key]=load(p,raw=True)[0];rows[key]={'owner':'world'};paths.append(p)
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis'])
    d=json.loads((OUT/'elbow_muscles/FJ1477/design.json').read_text());A,B=np.array(d['pad_centres_display_mm'])+placements['FJ1477'];r=B-E;trials=[];selected=None
    candidates=[]
    for phase in (np.arange(-44,-64,-2) if extended else np.arange(-42,-76,-2)):
        rotation=Rotation.from_rotvec(axis*np.radians(phase)).as_matrix()
        for axial in ([0,-8,8,-16,16,-20,20,-24,24,-28,28,-32,32] if extended else [0,-4,4,-8,8,-12,12,-16,16]):
            for scale in ([1.,1.1,1.2,1.3] if extended else [1.]):
                axialpart=axis*(r@axis);radial=r-axialpart
                delta=rotation@(radial*scale)+axialpart-r+axis*axial
                candidates.append((float(delta@delta),float(phase),axial,scale,delta))
    for _,phase,axial,scale,delta in sorted(candidates,key=lambda x:x[:4]):
        moments=[]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();P=E+R@(B+delta-E)
            moments.append(float(np.cross(P-E,unit(A-P))@axis))
        trial={'endpoint_phase_change_deg':phase,'axial_shift_mm':axial,'radial_scale':scale,'delta_display_mm':delta.tolist(),'delta_norm_mm':float(np.linalg.norm(delta)),'minimum_extension_chord_moment_magnitude_mm':float(min(-np.array(moments))),'first_wrong_chord_sign_deg':next((i for i,m in enumerate(moments) if m>=0),None)}
        if min(-np.array(moments))<2:
            trials.append(trial);continue
        moved={k:s.translate(delta) for k,s in moving.items()};first=None
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];hits=[]
            for key,s in moved.items():
                current=s.transform(T)
                for other,body in fixed.items():
                    target=body.transform(T) if rows[other]['owner']=='elbow' else body
                    v=intersection(current,target)
                    if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
            if hits:first={'angle_deg':angle,'hits':hits};break
        trial['first_hard_collision']=first;trials.append(trial)
        print('extensor endpoint',phase,axial,'first',first,'momentmargin',trial['minimum_extension_chord_moment_magnitude_mm'],flush=True)
        if first is None:selected=trial;break
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    report={'scope':'finite minimum-displacement endpoint family with negative chord moment screening and real distal hardware clearance, not compliant muscle validation','angle_sample_deg':1,'trials':trials,'selected':selected,'source_sha256':hashes,'old_yoke_reusable':False,'bone_centres_or_elbow_axis_changed':False,'actual_neutral_TPU_regenerated':False,'TPU_actuation_verified':False,'frozen':False}
    (dest/'endpoint_search.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('selected',selected,flush=True)

if __name__=='__main__':run()
