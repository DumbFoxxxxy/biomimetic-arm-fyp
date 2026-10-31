"""M2 real proximal grounding adapters; distal compliant output remains open.
Source anatomy, E1 skeleton and M1 standalone modules are preserved separately.
"""
from parameterize_elbow_muscles import *
from build_AC_connection import load
from build_SC_connection import rod,write,unit
from build_GH_connection import clean
import os
MUSCLES=OUT/'elbow_muscles'
SKELETON=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
DEST=MUSCLES/'mounted_elbow_candidate'

def line_bone_points(v,f,G,F,h,n):
    local=(v-G)@F;poly=section(local,f,h)[0];centre=pole(poly)[0];direction=F[:,:2].T@n;direction=direction/np.linalg.norm(direction);crossings=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        matrix=np.column_stack([direction,-(b-a)])
        if abs(np.linalg.det(matrix))<1e-10:continue
        t,s=np.linalg.solve(matrix,a-centre)
        if -.000001<=s<=1.000001:crossings.append(float(t))
    positive=[t for t in crossings if t>0];negative=[t for t in crossings if t<0]
    if not positive or not negative:raise RuntimeError('Bone mounting ray not bracketed')
    point=lambda t:G+F@np.r_[centre+direction*t,h]
    return point(min(positive)),point(max(negative))

def run():
    DEST.mkdir(exist_ok=True);definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'));G=np.array(definition['centres_source_mm']['GH'])-SHIFT
    hum,hv,hf=load(SKELETON/'humerus_GH_elbow_fork.obj');original_hum=hum
    manifest=json.loads((SKELETON/'whole_arm_STEP_check.json').read_text(encoding='utf8'));motion=json.loads((SKELETON/'bone_motion_audit.json').read_text(encoding='utf8'))
    instances=[];parts={};adapters={};placements=[];module_hard={};bone_hardware={};holes=[]
    for fid,distance,normal in [('FJ1486',16.,[0,-1,0]),('FJ1477',14.,[0,1,0])]:
        folder=MUSCLES/fid;data=json.loads((folder/'design.json').read_text(encoding='utf8'));F=np.array(data['basis_columns_display']);z=F[:,2];x=F[:,0];n=np.array(normal,dtype=float);n=unit(n-z*(n@z));offset=np.array(normal,dtype=float)*distance
        pad=np.array(data['pad_centres_display_mm'][0])+offset;cap_outer=pad-z*6.15;plate_centre=cap_outer-z*1.65
        station=float((plate_centre-G)@z);stations=[station-32,station-16] if station>90 else [52.,62.]
        flange=rod(cap_outer-z*3.15,cap_outer-z*.15,11.)
        if fid=='FJ1477':
            # The full round proximal flange contacts the scapula in axial
            # rotation. Trim only its unloaded posterior rim; both M3 washer
            # footprints and the boneward post attachment remain supported.
            perpendicular=unit(np.cross(z,n));frame=np.stack([perpendicular,z,n],1)
            keep=mf.Manifold.cube([200,200,200]).translate([-100,-100,-194.5]).transform(np.c_[frame,plate_centre])
            flange=flange^keep
            bevel_normal=n-z*(1.9/3);bevel_limit=4.55/np.linalg.norm(bevel_normal);bevel_normal=unit(bevel_normal);a=unit(np.cross(z,bevel_normal));b=np.cross(bevel_normal,a)
            bevel=mf.Manifold.cube([200,200,200]).translate([-100,-100,bevel_limit-200]).transform(np.c_[np.stack([a,b,bevel_normal],1),plate_centre]);flange=flange^bevel
        # Open the conduit-neck clearance toward the module exterior. Clamping
        # remains via the two retained M3 footprints, not via a closed neck ring.
        slot_frame=np.stack([unit(np.cross(z,n)),z,n],1)
        port_slot=mf.Manifold.cube([6.6,40,40]).translate([-3.3,-20,0]).transform(np.c_[slot_frame,plate_centre]);flange=flange-port_slot
        towards=unit((G+z*station-plate_centre)-z*((G+z*station-plate_centre)@z));root=plate_centre+towards*8;post=root-z*3.35;adapter=flange+rod(root,post,3.)
        anchors=[];brace_paths=[]
        for i,(h,side) in enumerate(zip(stations,[-1,1])):
            front,back=line_bone_points(hv,hf,G,F,h,n);qfront=front+n*3;qback=back-n*3
            hum=hum+rod(front-n*3,qfront,4.5)+rod(back+n*3,qback,4.5)
            hole=rod(qback-n*20,qfront+n*20,1.7);holes.append(hole)
            ear=rod(qfront+n*.15,qfront+n*3.15,5.)-hole;ear_centre=qfront+n*1.65
            side_start=post+x*side*25;side_end=ear_centre+x*side*25
            path=[post,side_start,side_end,ear_centre]
            adapter=adapter+ear+mf.Manifold.batch_boolean([rod(a,b,3.) for a,b in zip(path,path[1:])],mf.OpType.Add);brace_paths.append([p.tolist() for p in path])
            # Continuous head/washer approach recess, not only final-state seats.
            adapter=adapter-rod(qfront+n*3.0,qfront+n*20,3.75)
            head=qfront+n*3.65;required=float((head-qback)@n)+3.9
            length=next((L for L in [12,16,20,25,30,35,40,45,50] if L>=required),None)
            if length is None:raise RuntimeError('Mount bolt too long')
            name=fid+'_bone_'+str(i)
            bone_hardware[name+'_bolt']=rod(head,head-n*length,1.5)+rod(head,head+n*3,2.75)
            bone_hardware[name+'_front_washer']=rod(qfront+n*3.15,head,3.5)-hole
            bone_hardware[name+'_back_washer']=rod(qback,qback-n*.5,3.5)-hole
            xx=unit(np.cross(-n,z));yy=np.cross(-n,xx);frame=np.stack([xx,yy,-n],1)
            bone_hardware[name+'_nut']=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[frame,qback-n*.5])-rod(qback+n*4,qback-n*7,1.6)
            anchors.append({'station_below_GH_mm':h,'front_source_bone_mm':front.tolist(),'back_source_bone_mm':back.tolist(),'front_boss_face_mm':qfront.tolist(),'back_boss_face_mm':qback.tolist(),'bolt_axis':n.tolist(),'nominal_M3_bolt_length_mm':length})
        # Shared cap fasteners now also clamp the 3mm adapter flange.
        overrides={}
        for side,sign in [('left',-1),('right',1)]:
            q=pad+x*sign*7;head=q-z*9.8;name=fid+'_proximal_bolt_'+side
            head_height=1.65 if fid=='FJ1477' else 3.;head_radius=2.85 if fid=='FJ1477' else 2.75
            overrides['proximal_bolt_'+side]=name;parts[name]=rod(head,head+z*16,1.5)+rod(head,head-z*head_height,head_radius)
            source=load(folder/('proximal_head_washer_'+side+'.obj'))[0]
            name=fid+'_proximal_head_washer_'+side;overrides['proximal_head_washer_'+side]=name;parts[name]=source.translate(offset-z*3.15)
        flange_bores=[rod(pad+x*side*7-z*20,pad+x*side*7+z*8,1.7) for side in [-1,1]]
        adapter=adapter-rod(pad-z*25,pad+z*5,3.3)-mf.Manifold.batch_boolean(flange_bores+holes,mf.OpType.Add)
        adapters[fid+'_bone_adapter']=adapter
        for name in data['components']:
            key=fid+'_'+name
            if name in overrides:continue
            owner='GH' if name.startswith('proximal_') else 'soft_unvalidated' if name=='TPU_muscle' else 'distal_unvalidated'
            instances.append({'key':key,'file':str((folder/(name+'.obj')).relative_to(ROOT)),'translation_display_mm':offset.tolist(),'owner':owner,'kind':'TPU' if name=='TPU_muscle' else 'rigid','source_module':fid})
            if name!='TPU_muscle':module_hard[key]=load(folder/(name+'.obj'))[0].translate(offset)
        for name,s in parts.items():
            if name.startswith(fid+'_'):module_hard[name]=s
        placements.append({'id':fid,'translation_display_mm':offset.tolist(),'reason':'M1 rigid-clear screen plus mounting/fastener clearance; extensor moved further after initial soft/hardware interference. Provisional, no motion/TPU pass.','proximal_pad_centre_mm':pad.tolist(),'shared_cap_bolts':'M3x16 replaces standalone M3x12 at proximal end only','adapter_flange_thickness_mm':3.,'cap_flange_clearance_mm':.15,'bone_anchors':anchors,'brace_paths_display_mm':brace_paths})
    hum=hum-mf.Manifold.batch_boolean(holes,mf.OpType.Add);hum=clean(hum,'humerus_M2_mounts')
    # Explicit seats against neighboring hard geometry; collisions still checked later.
    mount_seats=[]
    for name,s in adapters.items():
        clearance=mf.Manifold.batch_boolean([hum]+[hum.translate(d*sign*.15) for d in np.eye(3) for sign in [-1,1]],mf.OpType.Add)
        s=s-clearance
        for other,solid in [*module_hard.items(),*bone_hardware.items()]:
            a=np.array(s.bounding_box());b=np.array(solid.bounding_box())
            if np.any(np.maximum(a[:3]-b[3:],b[:3]-a[3:])>=0):continue
            if (s^solid).volume()>.00001:
                before=float((s^solid).volume());cut=mf.Manifold.batch_boolean([solid]+[solid.translate(d*sign*.15) for d in np.eye(3) for sign in [-1,1]],mf.OpType.Add);s=s-cut;mount_seats.append({'adapter':name,'neighbor':other,'overlap_before_mm3':before,'feature':'local hard-part/head/washer clearance, not a collision exclusion'})
        adapters[name]=clean(s,name)
    parts={'humerus_M2_mounts':hum,**parts,**adapters,**bone_hardware}
    rows=[]
    for name,s in parts.items():
        s=clean(s,name);parts[name]=s;rows.append(write(DEST/(name+'.stl'),s))
        instances.append({'key':'FJ3368' if name=='humerus_M2_mounts' else name,'file':str((DEST/(name+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':'GH','kind':'rigid','generated':True})
    for key,record in manifest['parts'].items():
        if key=='FJ3368':continue
        path=ROOT/record['source_STEP'];obj=path.parent.parent/(record['part']+'.obj')
        instances.append({'key':key,'file':str(obj.relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':motion['moving_owners'].get(key,'world'),'kind':'rigid'})
    report={'scope':'M2 neutral proximal grounding candidate; distal tendon transmission, soft deformation and whole-arm drive unfinished','placements':placements,'instances':instances,'generated_components':{r['name']:r['components'] for r in rows},'humerus_delta_from_E1_mm3':{'removed':max(0,float((original_hum-hum).volume())),'added':max(0,float((hum-original_hum).volume()))},'material':'rigid bracket/bone materials unselected; nominal metal fastener envelopes, no thread/preload/strength validation','frozen':False,'pending':['all-pair neutral rigid and soft-contact classification','mount assembly and source-to-CAD consistency','moving proximal hardware vs torso/arm sweeps','distal guides, tendon routing/moment arms, required muscle stroke','actual TPU compression, printed fit, crimp grip, fatigue and loads']}
    report['local_mount_seats']=mount_seats
    report['anchor_revision']={'module':'FJ1477','first_station_below_GH_mm':{'before':50,'after':52},'reason':'2mm axial shift clears bone-bolt front washer approach past proximal module nut; joint centre unchanged'}
    report['adapter_port_slots']={'width_mm':6.6,'opens':'outward along each module mounting normal','purpose':'side insertion around6mm conduit neck, retain both shared M3 connections'}
    report['extensor_flange_revision']={'original_radius_mm':11,'posterior_extent_at_distal_face_mm':5.5,'posterior_extent_at_proximal_face_mm':3.6,'retained_M3_washer_radius_mm':3.5,'minimum_nominal_edge_beyond_washer_mm':.1,'reason':'full round flange first contacts scapula at~63deg; flat trim~80.6deg; preserve all washer support footprint and bevel only remaining unused posterior edge. Strength not validated.'}
    report['extensor_shared_screw_revision']={'family':'M3x16 button-head candidate, ISO7380-1 geometry envelope','maximum_head_diameter_mm':5.7,'maximum_head_height_mm':1.65,'modeled_as':'conservative full cylindrical head envelope; drive/thread not modeled','reference':'https://www.brightonbest.com/download/pds/PDS_Metric_Socket_Button_Head_Cap_Screws.pdf','reference_page':1,'reason':'trimmed flange gives~80.6deg positive axial prefix, next restriction is proximal screw head. No inherited load-capacity claim.'}
    (DEST/'assembly.json').write_text(json.dumps(report,indent=2),encoding='utf8');(DEST/'cad_input.json').write_text(json.dumps({'units':'mm','parts':rows}),encoding='utf8')
    print('M2 generated',len(parts),'instances',len(instances),'components',report['generated_components'],flush=True)

if __name__=='__main__':run()
