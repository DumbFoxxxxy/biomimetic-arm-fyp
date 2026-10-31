"""Direct TPU-end-to-forearm force path, no external passive steel loop.

One real material piece joins both distal caps to two existing ulna seats.
Caps keep their actual M3 clamps, wire bores and crimp seats. Muscle elasticity
is not prescribed by the rigid model and remains unverified.
"""
from build_forearm_receivers import ROOT,OUT,BASE,load,mf,np,rod,unit,write,intersection
import json,hashlib,itertools,sys
M2=OUT/'elbow_muscles/mounted_elbow_candidate'
SHIFTED='--shifted' in sys.argv
EXTENSOR_REANCHOR='--extensor-reanchor' in sys.argv
DEST=OUT/'direct_muscle_interfaces'/('mechanical_yoke_shifted' if SHIFTED else 'mechanical_yoke_zero')
if '--access-relief' in sys.argv:
    assert SHIFTED,'access relief belongs to the shifted candidate'
    DEST=OUT/'direct_muscle_interfaces/yoke_access_relief_study'
if EXTENSOR_REANCHOR:
    assert SHIFTED and '--access-relief' not in sys.argv
    DEST=OUT/'direct_muscle_interfaces/mechanical_yoke_reanchored_pair'

def muscle_folder(fid):
    if EXTENSOR_REANCHOR and fid=='FJ1477':return OUT/'direct_muscle_interfaces/extensor_reanchor_extended/muscle/FJ1477'
    if SHIFTED and fid=='FJ1486':return OUT/'direct_muscle_interfaces/shifted_muscle/FJ1486'
    return OUT/'elbow_muscles'/fid

def path(points,r=3):
    return mf.Manifold.batch_boolean([rod(a,b,r) for a,b in zip(points,points[1:])]+[mf.Manifold.sphere(r,16).translate(p) for p in points[1:-1]],mf.OpType.Add)

def run():
    assert '--access-relief' not in sys.argv,'Use materialize_direct_yoke_relief.py for the independent relief candidate'
    DEST.mkdir(exist_ok=True);a=json.loads((M2/'assembly.json').read_text());placements={p['id']:np.array(p['translation_display_mm']) for p in a['placements']};folders={fid:OUT/'elbow_muscles'/fid for fid in placements}
    folders={fid:muscle_folder(fid) for fid in placements}
    for r in a['instances']:
        for fid in placements:
            original=OUT/'elbow_muscles'/fid
            if (ROOT/r['file']).parent==original:r['file']=str((folders[fid]/(r['key'][7:]+'.obj')).relative_to(ROOT))
    parent4=json.loads((BASE/'assembly.json').read_text());bone_rows=[r.copy() for r in parent4['instances'] if r['key'].startswith(('bone_0_','bone_1_'))]
    for r in a['instances']:
        if r['key']=='FJ3391':r['file']=str((BASE/'forearm_spotfaced.obj').relative_to(ROOT));r['translation_display_mm']=[0,0,0]
    a['instances']+=bone_rows;rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};definition=json.loads((BASE/'design.json').read_text());n=np.array(definition['bone_anchors'][0]['bolt_axis']);down=unit(np.array(definition['zero_forearm_direction']));side=unit(np.cross(n,down));caps={};necks=[];cap_voids=[];roots=[];ears=[]
    for fid in placements:
        d=json.loads((folders[fid]/'design.json').read_text());F=np.array(d['basis_columns_display']);p=np.array(d['pad_centres_display_mm'][1])+placements[fid];x,y,z=F.T;caps[fid]=world[fid+'_distal_rigid_cap'];sign=1 if fid=='FJ1486' else -1;inner=p+F@np.array([0,sign*7.5,4.15]);outer=p+F@np.array([0,sign*14.5,4.15]);necks.append((fid,inner,outer));cap_voids+=[rod(p-z*16,p+z*16,.8),rod(p+z*2.95,p+z*6.25,1.4),rod(p+y*4+z*4.15-x*11,p+y*4+z*4.15+x*11,1.)]
        for direction in [-1,1]:q=p+x*direction*7;cap_voids.append(rod(q-z*12,q+z*12,1.7))
        if EXTENSOR_REANCHOR and fid=='FJ1486':
            q=p+x*7;cap_voids.append(rod(q+z*6.15,q+z*42,3.5,48))
    bone_voids=[];updated=[]
    for i,anchor in enumerate(definition['bone_anchors']):
        face=np.array(anchor['front_spotface_mm']);root=face+n*1.5;roots.append(root);ears.append(rod(face,face+n*3,5)-rod(face-n*8,face+n*20,1.7));bone_voids+=[rod(face-n*2,face+n*30,1.7),rod(face+n*3,face+n*30,3.75)]
        for suffix in ['bolt','front_washer']:
            key=f'bone_{i}_{suffix}';s=world[key].translate(-n*.15);updated.append(write(DEST/(key+'.stl'),s));world[key]=s;rows[key]['file']=str((DEST/(key+'.obj')).relative_to(ROOT));rows[key]['translation_display_mm']=[0,0,0]
    against={k:s for k,s in world.items() if not k.endswith('_distal_rigid_cap')};trials=[];chosen=None
    # Check bony context as well as all manufactured neutral hard/TPU geometry.
    refs=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in refs['models']:
        if r['kind']=='bone_reference' and r['id'] not in against:against['reference_'+r['id']]=load(OUT/'repaired_supports'/(r['id']+'.obj'),raw=True)[0]
    voids=mf.Manifold.batch_boolean(cap_voids+bone_voids,mf.OpType.Add);original_forearm=world['FJ3391'];seat_region=mf.Manifold.batch_boolean([mf.Manifold.sphere(14,24).translate(root+down*4) for root in roots],mf.OpType.Add)
    for reach,dodge,extra in itertools.product([10,16,22,28,36],[-18,18],[0,10,20]):
        exterior=[root+down*8+n*reach+side*dodge for root in roots];centre=np.mean(exterior,axis=0);pieces=ears+list(caps.values());paths=[]
        for root,end in zip(roots,exterior):
            points=[root,root+down*8,root+down*8+n*reach,end];pieces.append(path(points));paths.append([p.tolist() for p in points])
        pieces.append(rod(exterior[0],exterior[1],3))
        for fid,inner,outer in necks:
            outward=-1 if fid=='FJ1486' else 1;stage=centre.copy();stage[1]=outer[1]+outward*extra;stage[2]=outer[2]-10;points=[centre,stage,outer];pieces+=[path(points),rod(inner,outer,2),mf.Manifold.sphere(3,16).translate(outer)];paths.append([p.tolist() for p in points])
        yoke=(mf.Manifold.batch_boolean(pieces,mf.OpType.Add)-voids).as_original().simplify(.001);components=len([s for s in yoke.decompose() if s.volume()>.001]);hits=[]
        removed=original_forearm^yoke;outside=max(0,float((removed-seat_region).volume()));fitted=(original_forearm-yoke).as_original().simplify(.001);fitted_components=len([s for s in fitted.decompose() if s.volume()>.001]);trial_against={**against,'FJ3391':fitted}
        for k,s in trial_against.items():
            v=intersection(yoke,s)
            if v>.001:hits.append({'part':k,'volume_mm3':v,'kind':'TPU_review' if k in rows and rows[k]['kind']=='TPU' else 'rigid'})
        trial={'reach_mm':reach,'side_dodge_mm':dodge,'extra_clearance_mm':extra,'components':components,'neutral_interferences':hits,'paths_display_mm':paths,'controlled_forearm_seat_removed_mm3':float(removed.volume()),'removed_outside_two_14mm_mount_regions_mm3':outside,'forearm_positive_components':fitted_components};trials.append(trial)
        if components==1 and not hits and outside<.001 and removed.volume()<100 and fitted_components==1:chosen=(yoke,trial,fitted);break
    report={'scope':'one-piece direct muscle cap/forearm yoke, neutral manufactured geometry only, not soft actuation or installed motion acceptance','shifted_flexor':SHIFTED,'trials':trials,'selected':chosen[1] if chosen else None,'bones_and_joint_centres_changed':False,'original_cap_ports_preserved':True,'material_and_strength_verified':False,'frozen':False};(DEST/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('direct yoke neutral candidates',len(trials),'selected',bool(chosen),flush=True)
    if chosen is None:return
    yoke,trial,fitted=chosen;updated.append(write(DEST/'direct_muscle_yoke.stl',yoke));updated.append(write(DEST/'direct_forearm_seated.stl',fitted));rows['FJ3391']['file']=str((DEST/'direct_forearm_seated.obj').relative_to(ROOT));a['instances']=[r for r in a['instances'] if not r['key'].endswith('_distal_rigid_cap')];a['instances'].append({'key':'direct_muscle_yoke','file':str((DEST/'direct_muscle_yoke.obj').relative_to(ROOT)),'translation_display_mm':[0,0,0],'kind':'rigid','owner':'elbow'})
    uncrimped=[]
    for r in a['instances']:
        if any(r['key'].startswith(fid+'_distal_') or r['key']==fid+'_active_wire_crimp_barrel' for fid in placements):r['owner']='elbow'
        if r['key'].endswith('_active_wire_crimp_barrel'):r['nominal_packing_only']=True;r['nominal_pose_owner']='elbow';r['owner']='distal_unvalidated';uncrimped.append(r['key'])
    a.update(scope='direct TPU bone-interface candidate, actual rigid cap/yoke constraint; TPU elastic deformation unverified',frozen=False,uncrimped_nominal_stops=uncrimped);(DEST/'assembly.json').write_text(json.dumps(a,indent=2),encoding='utf8');(DEST/'cad_input.json').write_text(json.dumps({'units':'mm','parts':updated}),encoding='utf8');print('direct candidate manufactured instances',len(a['instances']),flush=True)

if __name__=='__main__':run()
