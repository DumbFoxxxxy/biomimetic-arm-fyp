"""Finite endpoint relocation search, with real distal hardware and bones.

Positions are candidate changes to TPU end geometry, NOT bone joint shifts.
The endpoint chord moment is only a screening approximation, not a muscle law.
"""
from probe_direct_end_hardware import *
from build_forearm_receivers import unit
import itertools

def run():
    a=json.loads((M2/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);down=unit(np.array(j['zero_forearm_direction']));side=np.cross(axis,down);all_distal={k:s for k,s in world.items() if any(k.startswith(fid+'_distal_') or k==fid+'_active_wire_crimp_barrel' for fid in ['FJ1486','FJ1477'])};flex={k:s for k,s in all_distal.items() if k.startswith('FJ1486_')};ext={k:s for k,s in all_distal.items() if k.startswith('FJ1477_')};others={k:s for k,s in world.items() if k not in all_distal and rows[k]['kind']!='TPU'}
    sourcepaths=[M2/'assembly.json',BASE/'design.json']+[ROOT/r['file'] for r in rows.values()]
    for r in json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))['models']:
        if r['kind']=='bone_reference' and r['id'] not in others:
            p=OUT/'repaired_supports'/(r['id']+'.obj');key='reference_'+r['id'];others[key]=load(p,raw=True)[0];rows[key]={'owner':'world'};sourcepaths.append(p)
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcepaths};d=json.loads((OUT/'elbow_muscles/FJ1486/design.json').read_text());offset=np.array(next(p for p in a['placements'] if p['id']=='FJ1486')['translation_display_mm']);A,B=np.array(d['pad_centres_display_mm'])+offset
    trials=[];selected=None
    parameters=sorted(itertools.product([0,-4,4,-8,8,-12,12,-16,16,-20,20,-24,24,-28,28,-32,32,-40,40],[0,4,8,12,16]),key=lambda p:(p[0]**2+p[1]**2,p))
    for lateral,posterior in parameters:
        delta=axis*lateral-side*posterior;moments=[]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();P=E+R@(B+delta-E);moments.append(float(np.cross(P-E,unit(A-P))@axis))
        row={'lateral_offset_mm':lateral,'posterior_offset_mm':posterior,'delta_display_mm':delta.tolist(),'offset_norm_mm':float(np.linalg.norm(delta)),'minimum_chord_force_moment_arm_mm':min(moments),'first_wrong_chord_sign_sample_deg':next((i for i,m in enumerate(moments) if m<=0),None)}
        if row['first_wrong_chord_sign_sample_deg'] is not None:trials.append(row);continue
        moved={k:s.translate(delta) for k,s in flex.items()};cross=[]
        for k,s in moved.items():
            for b,body in ext.items():
                v=intersection(s,body)
                if v>.001:cross.append({'pair':[k,b],'volume_mm3':v})
        first=None
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];hits=[]
            for k,s in moved.items():
                candidate=s.transform(T)
                for b,body in others.items():
                    other=body.transform(T) if rows[b]['owner']=='elbow' else body;v=intersection(candidate,other)
                    if v>.001:hits.append({'pair':[k,b],'volume_mm3':v})
            if hits:first={'elbow_angle_deg':angle,'hits':hits};break
        row['first_hard_collision']=first;row['cross_module_neutral_collisions']=cross;trials.append(row);print('direct flex endpoint',lateral,posterior,'first collision',first,flush=True)
        if first is None and not cross:selected=row;break
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    report={'scope':'minimum sampled offset in this finite family passing real hard-hardware150deg and positive chord-force screening; not validated TPU deformation or direct adapter assembly','trials':trials,'selected':selected,'source_sha256':hashes,'step_deg':1,'bones_or_joint_centres_changed':False,'actual_muscle_geometry_regenerated':False,'drive_verified':False,'frozen':False};(DEST/'flexor_attachment_search.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('selected candidate',selected,flush=True)

if __name__=='__main__':run()
