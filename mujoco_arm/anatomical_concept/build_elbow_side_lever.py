"""M4 removable forearm-side cable lever, independent from rejected M3 rotor.

Actual ulna-referenced mounting seats, through bolts and tendon eyelet. This
does not invent a fixed guide or establish muscle shortening/deformation.
"""
from build_elbow_muscle_mounts import ROOT,OUT,DEST as MOUNT,load,rod,unit,mf,np
from build_SC_connection import write
from parameterize_elbow_muscles import section,pole
from build_reference import SHIFT,read_obj
import json,hashlib,sys

DEST=MOUNT/'side_lever_candidate'

def run():
    DEST.mkdir(exist_ok=True)
    skeleton=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
    elbow=json.loads((skeleton/'assembly.json').read_text(encoding='utf8'))
    E=np.array(elbow['centre_source_mm'])-SHIFT;axis=unit(np.array(elbow['axis_source']))
    bridge=np.array(elbow['forearm_bridges_display_mm'][0]);down=unit(bridge.mean(0)-E);down=unit(down-axis*(down@axis));y=np.cross(down,axis);F=np.stack([axis,y,down],1)
    inventory=json.loads((OUT/'inventory.json').read_text(encoding='utf8'))
    source=ROOT/next(r['source'] for r in inventory['models'] if r['id']=='FJ3391')
    uv,uf=read_obj(str(source));ulocal=(uv-SHIFT-E)@F
    forearm,v,f=load(skeleton/'forearm_elbow_hub.obj');local=(v-E)@F;stations=[]
    for h in [42.,58.]:
        centre=pole(section(ulocal,uf,h)[0])[0];crossings=[]
        for poly in section(local,f,h):
            for a,b in zip(poly,np.roll(poly,-1,axis=0)):
                if abs(b[1]-a[1])<1e-10:continue
                t=(centre[1]-a[1])/(b[1]-a[1])
                if -1e-8<=t<=1+1e-8:crossings.append(float(a[0]+t*(b[0]-a[0])))
        lower=max(x for x in crossings if x<centre[0]);upper=min(x for x in crossings if x>centre[0])
        stations.append({'station_mm':h,'ulnar_section_pole_local_mm':centre.tolist(),'axis_intersections_local_mm':[lower,upper]})
    variants=[]
    phase_search=json.loads((DEST/'eyelet_phase_search.json').read_text(encoding='utf8')) if '--phase-search' in sys.argv else None
    for sign in [-1,1]:
        n=axis*sign;anchors=[];cuts=[];hardware={};ears=[];roots=[]
        for i,s in enumerate(stations):
            h=s['station_mm'];lower,upper=s['axis_intersections_local_mm'];cy=s['ulnar_section_pole_local_mm'][1]
            front=E+F@np.array([upper if sign>0 else lower,cy,h]);back=E+F@np.array([lower if sign>0 else upper,cy,h])
            front_face=front-n*2;back_face=back+n*2
            # Recessed annular spotfaces remove material; no integral bosses
            # grow onto the forearm before it is seated in the elbow.
            hole=rod(back-n*30,front+n*30,1.7)
            cuts+=[hole,rod(front_face,front+n*20,5.15),rod(back_face,back-n*20,3.65)]
            ear=rod(front_face+n*.15,front_face+n*3.15,5)-hole;ears.append(ear)
            root=front_face+n*1.65;roots.append(root)
            head=front_face+n*3.65;required=float((head-back_face)@n)+3.9
            length=next(L for L in [12,16,20,25,30,35,40,45,50] if L>=required)
            name='bone_'+str(i)
            hardware[name+'_bolt']=rod(head,head-n*length,1.5)+rod(head,head+n*3,2.75)
            hardware[name+'_front_washer']=rod(front_face+n*3.15,head,3.5)-hole
            hardware[name+'_back_washer']=rod(back_face,back_face-n*.5,3.5)-hole
            frame=np.stack([down,np.cross(-n,down),-n],1)
            hardware[name+'_nut']=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[frame,back_face-n*.5])-rod(back_face+n*4,back_face-n*7,1.6)
            anchors.append({'station_mm':h,'front_surface_mm':front.tolist(),'back_surface_mm':back.tolist(),'front_spotface_mm':front_face.tolist(),'back_spotface_mm':back_face.tolist(),'bolt_axis':n.tolist(),'nominal_M3_length_mm':length})
        revised=(forearm-mf.Manifold.batch_boolean(cuts,mf.OpType.Add)).as_original().simplify(.01)
        assert revised.status()==mf.Error.NoError
        for radius in [10.,15.,20.]:
            name=f'axis_{"positive" if sign>0 else "negative"}_r{int(radius)}';folder=DEST/name;folder.mkdir(exist_ok=True)
            end=E+axis*sign*56+down*radius;opposite=E+axis*sign*56-down*radius;spine_end=E+axis*sign*56+down*58
            phases=None
            if phase_search:
                found=next(r for r in phase_search['variants'] if r['variant']==name)['muscles']
                if all(r['selected_diagnostic_candidate'] for r in found):
                    phases={r['muscle']:r['selected_diagnostic_candidate'] for r in found}
                    end=np.array(phases['FJ1486']['eyelet_display_mm']);opposite=np.array(phases['FJ1477']['eyelet_display_mm'])
            spine_base=E+axis*sign*56+down*radius
            lever=mf.Manifold.batch_boolean(ears+[rod(spine_base,spine_end,3.),rod(opposite,spine_base,3.),rod(end,spine_base,3.)] if np.linalg.norm(end-spine_base)>1e-6 else ears+[rod(opposite,spine_end,3.)],mf.OpType.Add)
            brace_paths=[]
            for root in roots:
                turn=root+down*8
                outer=turn+n*float((E+axis*sign*56-turn)@n)
                path=[root,turn,outer,end]
                lever=lever+mf.Manifold.batch_boolean([rod(a,b,3.) for a,b in zip(path,path[1:])],mf.OpType.Add);brace_paths.append([p.tolist() for p in path])
            # Eyelet is actual geometry; tendon knot/crimp and contact remain
            # unvalidated. The hole centre is only a routing reference.
            lever=lever+rod(end-axis*2,end+axis*2,5.)
            lever=lever-rod(end-axis*6,end+axis*6,.8)
            lever=lever+rod(opposite-axis*2,opposite+axis*2,5.)
            lever=lever-rod(opposite-axis*6,opposite+axis*6,.8)
            tendon_stops={}
            for label,point in [('flexor',end),('extensor',opposite)]:
                far=axis*sign;near=-far
                lever=lever-rod(point+far*2,point+far*8,1.6)-rod(point+near*1.8,point+near*8,2.)
                tendon_stops[label+'_tendon_stop']=rod(point+far*2.05,point+far*5.05,1.4)-rod(point+far*1.5,point+far*6,.55)
            for a in anchors:
                p=np.array(a['front_spotface_mm']);lever=lever-rod(p-n*2,p+n*20,1.7)-rod(p+n*3.0,p+n*20,3.75)
            lever=lever.as_original().simplify(.01)
            before_seat=max(0,float((revised^lever).volume()))
            seat=mf.Manifold.batch_boolean([lever]+[lever.translate(d*sign*.15) for d in np.eye(3) for sign in [-1,1]],mf.OpType.Add)
            fitted=(revised-seat).as_original().simplify(.01)
            assert len([c for c in fitted.decompose() if c.volume()>.001])==1,'contour seat disconnects forearm'
            parts={'forearm_spotfaced':fitted,'side_lever':lever,**hardware,**tendon_stops};rows=[]
            for key,solid in parts.items():
                assert solid.status()==mf.Error.NoError
                rows.append(write(folder/(key+'.stl'),solid))
            report={'stage':'M4 removable forearm side lever','variant':name,'axis_side':sign,'nominal_tendon_reference_radius_mm':radius,'eyelet_axis_station_mm':sign*56,'eyelet_centre_display_mm':end.tolist(),'eyelet_bore_radius_mm':.8,'nominal_wire_radius_mm':.3,'joint_centre_display_mm':E.tolist(),'joint_axis':axis.tolist(),'zero_forearm_direction':down.tolist(),'bone_anchors':anchors,'brace_paths':brace_paths,'parts':{r['name']:r['components'] for r in rows},'forearm_removed_mm3':max(0,float((forearm-revised).volume())),'no_new_integral_forearm_bosses':True,'source_forearm_sha256':hashlib.sha256((skeleton/'forearm_elbow_hub.obj').read_bytes()).hexdigest(),'force_path':'passive tendon loop at actual lever eyelet -> removable lever -> two M3 clamped spotfaces -> source ulna/radius rigid body; thread/preload/torque capacity not verified','muscle_or_guide_constraint_imposed':False,'drive_pass':False,'frozen':False}
            report.update(forearm_removed_mm3=max(0,float((forearm-fitted).volume())),contour_seat={'actual_brace_bone_overlap_before_mm3':before_seat,'additional_material_removed_mm3':max(0,float((revised-fitted).volume())),'construction':'actual lever and six +/-0.15mm Cartesian offsets subtracted locally; no collision pair excluded','positive_material_components':1,'strength_verified':False})
            report['eyelets_display_mm']={'distal':end.tolist(),'proximal':opposite.tolist()}
            if phases:report['eyelets_display_mm']={'flexor':end.tolist(),'extensor':opposite.tolist()}
            report['eyelet_phase_diagnostic']=phases
            report['free_span_entry_offset_towards_cap_mm']=2.4
            report['tendon_stops']={'crimp_barrel_outer_diameter_mm':2.8,'uncrimped_bore_diameter_mm':1.1,'length_mm':3.,'nominal_stop_to_back_face_gap_mm':.05,'eyelet_through_diameter_mm':1.6,'actual_crimp_grip_verified':False,'steel_entry_bend_verified':False}
            (folder/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');(folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':rows}),encoding='utf8');variants.append(report)
            print(name,'components',report['parts'],'removed',report['forearm_removed_mm3'],flush=True)
    (DEST/'variants.json').write_text(json.dumps({'scope':'manufactured rigid lever candidates; not solved tendon/muscle motion','variants':[r['variant'] for r in variants]},indent=2),encoding='utf8')
if __name__=='__main__':run()
