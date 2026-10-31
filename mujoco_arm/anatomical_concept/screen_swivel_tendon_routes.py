"""Real M5 ball centres -> real M4 rear-stop centres, straight route screen.

Muscle/cap centres are held at neutral only to screen geometry, NOT because a
rail or prescribed axial muscle law exists. No unsupported outlet kink is used.
"""
from build_distal_swivels import *
from scipy.spatial.transform import Rotation
def align(a,b):
    cross=np.cross(a,b);n=np.linalg.norm(cross);dot=float(np.clip(a@b,-1,1))
    if n<1e-10:return np.eye(3)
    return Rotation.from_rotvec(cross/n*np.arccos(dot)).as_matrix()
def run():
    assembly=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));rows={r['key']:r for r in assembly['instances']};source={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};lever=json.loads((BASE/'design.json').read_text(encoding='utf8'));E=np.array(lever['joint_centre_display_mm']);axis=np.array(lever['joint_axis']);modules={fid:json.loads((DEST/fid/'design.json').read_text(encoding='utf8')) for fid in ['FJ1486','FJ1477']}
    support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'));context=[]
    for r in support['models']:
        if r['kind']!='bone_reference' or r['id']=='FJ3290':continue
        b=load(OUT/'repaired_supports'/(r['id']+'.obj'),True)[0]
        for key,s in source.items():
            if '_swivel_' not in key:continue
            volume=intersection(s,b)
            if volume>.001:context.append({'pair':[key,r['id']],'volume_mm3':volume})
    samples=[];lengths={fid:[] for fid in modules}
    for angle in range(151):
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];world={k:s.transform(T) if rows[k]['owner']=='elbow' else s for k,s in source.items()};predicted={}
        for fid,eyelet in [('FJ1486','flexor'),('FJ1477','extensor')]:
            p=modules[fid]['selected'];C=np.array(p['ball_centre_display_mm']);u=np.array(p['neutral_terminal_axis']);anchor0=np.array(lever['eyelets_display_mm'][eyelet]);B=E+R@(anchor0-E)+axis*lever['axis_side']*5.5;direction=unit(B-C);rotation=align(u,direction)
            world[fid+'_swivel_ball_crimp_blank']=source[fid+'_swivel_ball_crimp_blank'].transform(np.c_[rotation,C-rotation@C]);wire=rod(C+direction*3.5,B,.3);lengths[fid].append(float(np.linalg.norm(B-C)))
            predicted[fid]={'C':C,'B':B,'wire':wire,'direction_change_deg':float(np.degrees(np.arccos(np.clip(u@direction,-1,1))))}
        row={'angle_deg':angle,'muscles':[]}
        for fid,p in predicted.items():
            metal_key=fid+'_swivel_ball_crimp_blank';wire_hits=[];terminal_hits=[];soft=[]
            for key,s in world.items():
                volume=intersection(p['wire'],s)
                if volume>.001:(soft if rows[key]['kind']=='TPU' else wire_hits).append({'part':key,'volume_mm3':volume})
                if key!=metal_key and rows[key]['kind']!='TPU':
                    volume=intersection(world[metal_key],s)
                    if volume>.001:terminal_hits.append({'part':key,'volume_mm3':volume})
            row['muscles'].append({'muscle':fid,'terminal_direction_change_deg':p['direction_change_deg'],'ball_centre_display_mm':p['C'].tolist(),'actual_stop_reference_display_mm':p['B'].tolist(),'straight_wire_penetrations':wire_hits,'terminal_metal_interferences':terminal_hits,'static_TPU_contact_flags':soft})
        samples.append(row)
        if angle%30==0:print('M5 actual-endpoint route screen',angle,[(r['muscle'],len(r['straight_wire_penetrations'])) for r in row['muscles']],flush=True)
    summary=[]
    for fid in modules:
        cases=[(s['angle_deg'],next(r for r in s['muscles'] if r['muscle']==fid)) for s in samples];desired=1 if fid=='FJ1486' else -1;moment=-np.gradient(lengths[fid],np.radians(np.arange(151)));wrong=np.flatnonzero(moment*desired<=0);hits=[{'angle_deg':angle,'hits':r['straight_wire_penetrations']} for angle,r in cases if r['straight_wire_penetrations']]
        summary.append({'muscle':fid,'held_neutral_cap':True,'maximum_required_terminal_rotation_deg':max(r['terminal_direction_change_deg'] for _,r in cases),'straight_wire_penetration_samples':len(hits),'first_wire_penetration':hits[0] if hits else None,'terminal_metal_interference_samples':sum(bool(r['terminal_metal_interferences']) for _,r in cases),'direct_endpoint_distance_mm':lengths[fid],'geometric_distance_excursion_mm':float(np.ptp(lengths[fid])),'minus_dL_dtheta_mm':moment.tolist(),'first_wrong_or_zero_direct_torque_sign_deg':int(wrong[0]) if len(wrong) else None})
    report={'scope':'actual spherical attachment and actual rear-stop references, no unsupported outlet-point bends. Cap centres held neutral for geometric screening only; not a muscle constraint/deformation or tension solution. Stops remain uncrimped, linear wire penetrations require real routing.','angular_step_deg':1,'rigid_threshold_mm3':.001,'new_terminal_vs_fixed_bony_context':context,'muscles':summary,'samples':samples,'assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest(),'drive_pass':False,'frozen':False}
    (DEST/'straight_route_screen.json').write_text(json.dumps(report,indent=2),encoding='utf8');print([{k:v for k,v in r.items() if k not in ['direct_endpoint_distance_mm','minus_dL_dtheta_mm']} for r in summary],flush=True)
if __name__=='__main__':run()
