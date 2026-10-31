"""Local rigid fitting and explicit soft-contact flags for direct end mounting.

This does not infer complete assembly from a rigid approximation to the TPU.
"""
from build_direct_muscle_yoke import *
ACCESS_STUDY='--access-study' in sys.argv
PRELOADED='--preloaded-sequence' in sys.argv
RELIEF_STUDY='--relief-study' in sys.argv

def run():
    a=json.loads((DEST/'assembly.json').read_text());rows={r['key']:r for r in a['instances']}
    if RELIEF_STUDY:
        rows['direct_muscle_yoke']['file']=str((OUT/'direct_muscle_interfaces/yoke_access_relief_study/direct_muscle_yoke.obj').relative_to(ROOT))
    world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};results=[];inputs=[DEST/'assembly.json']+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    def probe(label,moving,fixed,direction,length):
        if ACCESS_STUDY:
            if label not in ['FJ1486_proximal_left','FJ1486_distal_right']:return
            def check_offsets(offsets):
                soft=[]
                for offset in offsets:
                    for key in moving:
                        moved=world[key].translate(offset)
                        for other in fixed:
                            v=intersection(moved,world[other])
                            if v>.001:
                                hit={'offset_mm':offset.tolist(),'pair':[key,other],'volume_mm3':v}
                                if 'TPU' in [rows[key]['kind'],rows[other]['kind']]:soft.append(hit)
                                else:return {'rigid_clear':False,'first_rigid_limit':hit,'soft_flags':soft}
                return {'rigid_clear':True,'first_rigid_limit':None,'soft_flags':soft}
            trials=[];selected=None
            for stage in [3,5,8,10,12,15,18,20]:
                closing=[direction*d for d in np.linspace(stage,0,int(stage*10)+1)]
                closure=check_offsets(closing)
                if not closure['rigid_clear']:
                    trials.append({'axial_stage_mm':stage,'closure':closure});continue
                for reference in np.eye(3):
                    lateral=reference-direction*(reference@direction)
                    if np.linalg.norm(lateral)<.1:continue
                    lateral=unit(lateral)
                    for sign in [-1,1]:
                        side=lateral*sign
                        offsets=[side*d+direction*stage for d in np.r_[np.arange(40,5,-1),np.linspace(5,0,51)]]
                        result=check_offsets(offsets)
                        record={'axial_stage_mm':stage,'side_direction':side.tolist(),'side_entry_mm':40,'side':result,'closure':closure}
                        trials.append(record)
                        if result['rigid_clear']:
                            selected=record;break
                    if selected:break
                if selected:break
            results.append({'label':label,'moving':moving,'fixed':fixed,'trials':trials,'selected':selected})
            print(label,'staged access found',bool(selected),flush=True)
            return
        first=None;soft=[]
        for distance in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
            h=[]
            for key in moving:
                moved=world[key].translate(direction*distance)
                for other in fixed:
                    v=intersection(moved,world[other])
                    if v>.001:
                        hit={'pair':[key,other],'volume_mm3':v}
                        if 'TPU' in [rows[key]['kind'],rows[other]['kind']]:soft.append({'distance_mm':float(distance),**hit})
                        else:h.append(hit)
            if h and first is None:first={'distance_mm':float(distance),'hits':h}
        r={'label':label,'moving':moving,'direction':direction.tolist(),'rigid_clear':first is None,'first_rigid_limit':first,'soft_flags':soft,'final5mm_step_mm':.1,'other_step_mm':1};results.append(r);print(label,'rigid',first is None,'soft',len(soft),flush=True)
    j=json.loads((BASE/'design.json').read_text());n=np.array(j['bone_anchors'][0]['bolt_axis']);newbonehw=[f'bone_{i}_{s}' for i in [0,1] for s in ['bolt','front_washer','back_washer','nut']];muscle=lambda k:k.startswith(('FJ1486_','FJ1477_')) and '_bone_' not in k
    late=[f'FJ1486_bone_{i}_{suffix}' for i in [0,1] for suffix in ['bolt','front_washer']] if RELIEF_STUDY or EXTENSOR_REANCHOR else []
    fixed=[k for k in world if not muscle(k) and k not in newbonehw+['direct_muscle_yoke']+late]
    early_distal=['FJ1486_distal_bolt_right','FJ1486_distal_head_washer_right'] if PRELOADED else []
    if PRELOADED:
        d=json.loads((OUT/'direct_muscle_interfaces/shifted_muscle/FJ1486/design.json').read_text());z=np.array(d['basis_columns_display'])[:,2]
        probe('FJ1486_distal_right_into_bare_yoke',early_distal,['direct_muscle_yoke'],z,25)
    probe('direct_yoke_to_forearm_before_muscles',['direct_muscle_yoke']+early_distal,fixed,n,60);fixed+=['direct_muscle_yoke']+early_distal
    for i in [0,1]:
        prefix=f'bone_{i}';keys=[prefix+'_bolt',prefix+'_front_washer'];probe(prefix+'_front',keys,fixed,n,45);fixed+=keys;keys=[prefix+'_nut',prefix+'_back_washer'];probe(prefix+'_rear',keys,fixed,-n,35);fixed+=keys
    for fid,normal in [('FJ1486',[0,-1,0]),('FJ1477',[0,1,0])]:
        folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());z=np.array(d['basis_columns_display'])[:,2];normal=np.array(normal,dtype=float);normal=unit(normal-z*(normal@z));bolts=[fid+'_'+end+'_'+suffix+'_'+side for end in ['proximal','distal'] for suffix in ['bolt','head_washer'] for side in ['left','right']];crimp=fid+'_active_wire_crimp_barrel';module=[k for k in world if k.startswith(fid+'_') and '_bone_' not in k and k not in bolts+[crimp]];probe(fid+'_neutral_shape_side_entry_to_two_caps',module,fixed,normal,25);fixed+=module
        for end,sign in [('proximal',-1),('distal',1)]:
            for side in ['left','right']:
                keys=[fid+'_'+end+'_bolt_'+side,fid+'_'+end+'_head_washer_'+side]
                if PRELOADED and keys==early_distal:continue
                probe(fid+'_'+end+'_'+side,keys,fixed,z*sign,25);fixed+=keys
        probe(fid+'_crimp_stop_before_gripping_wire',[crimp],fixed,z,12);fixed.append(crimp)
    if RELIEF_STUDY or EXTENSOR_REANCHOR:
        placement=next(p for p in a['placements'] if p['id']=='FJ1486');axis=np.array(placement['bone_anchors'][0]['bolt_axis'])
        for i in [0,1]:
            keys=[f'FJ1486_bone_{i}_{s}' for s in ['bolt','front_washer']]
            probe(f'FJ1486_bone_{i}_front_after_caps',keys,fixed,axis,30);fixed+=keys
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    if ACCESS_STUDY:
        report={'scope':'two previously obstructed bolt/head-washer pairs, lateral staging then axial closure; same installed fixed-part sets; samples <=1mm and final5mm <=0.1mm; no tool-access or soft-deformation proof','paths':results,'source_sha256':hashes,'complete_arm_assembly_verified':False,'frozen':False}
        (DEST/'bolt_access_study.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    else:
        target=OUT/'direct_muscle_interfaces/yoke_access_relief_study' if RELIEF_STUDY else DEST
        (target/('preloaded_assembly_checks.json' if PRELOADED else 'assembly_checks.json')).write_text(json.dumps({'scope':'sampled direct prototype fitting; body neutral form held for insertion only, positive soft flags require compliant/physical checking, no automatic failure or pass','paths':results,'rigid_clear_count':sum(r['rigid_clear'] for r in results),'soft_flag_count':sum(len(r['soft_flags']) for r in results),'source_sha256':hashes,'complete_arm_assembly_verified':False,'frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
