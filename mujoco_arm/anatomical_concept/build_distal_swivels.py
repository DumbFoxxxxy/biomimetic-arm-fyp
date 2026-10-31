"""M5 actual split spherical tendon terminals on anatomical TPU end caps.

Ball/crimp blank and two-bolt retainer give a physical force path; material,
crimp grip, strength, compliant muscle motion and cable routes stay unverified.
"""
from build_elbow_side_lever import ROOT,OUT,MOUNT,load,rod,unit,mf,np,write
import json,hashlib,itertools,sys

BASE=MOUNT/'side_lever_candidate/axis_negative_r15'
DISTAL_LOAD='--distal-load' in sys.argv
DEST=BASE/('swivel_terminal_distal_load' if DISTAL_LOAD else 'swivel_terminal_candidate')

def intersection(a,b):
    aa=np.array(a.bounding_box());bb=np.array(b.bounding_box())
    if np.any(np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])>=0):return 0.
    return max(0.,float((a^b).volume()))

def stepped_bolt(head,u,v,w):
    # One continuous24-sided profile avoids microscopic coplanar Boolean seams
    # between a shaft and its head, which can invalidate a sewn BRep.
    vertices=[]
    for height,radius in [(-10,1.),(0,1.),(0,1.9),(1.6,1.9)]:
        vertices += [(head+u*height+radius*(v*np.cos(a)+w*np.sin(a))).tolist() for a in np.arange(24)*np.pi/12]
    vertices += [(head-u*10).tolist(),(head+u*1.6).tolist()];faces=[]
    for ring in range(3):
        for j in range(24):
            a=ring*24+j;b=ring*24+(j+1)%24;c=(ring+1)*24+(j+1)%24;d=(ring+1)*24+j;faces += [[a,b,c],[a,c,d]]
    for j in range(24):faces += [[(j+1)%24,j,96],[72+j,72+(j+1)%24,97]]
    result=mf.Manifold(mf.Mesh64(np.array(vertices,dtype=float),np.array(faces,dtype=np.uint64)));assert result.status()==mf.Error.NoError
    return result

def construct(cap,C,u,v,root):
    w=np.cross(u,v);frame=np.stack([v,w,u],1)
    def half(lo,hi):return mf.Manifold.cube([100,100,hi-lo]).translate([-50,-50,lo]).transform(np.c_[frame,C])
    sphere=mf.Manifold.sphere(5.2,32).translate(C);cavity=mf.Manifold.sphere(3.7,32).translate(C)
    # Explicit conical window: radius1.45 at the centre,25deg half angle.
    t=np.tan(np.radians(25));window=mf.Manifold.cylinder(12,1.45-t,1.45+11*t,48).transform(np.c_[frame,C-u])
    shell=sphere-cavity-window
    rear=cap+(shell^half(-20,-.15))+rod(root,C-u*4.5,2.5)
    lid=shell^half(.15,20);hardware={};eyes=[]
    for label,sign in [('left',-1),('right',1)]:
        p=C+v*sign*6.8;eyes.append(p)
        rear=rear+rod(p-u*3.5,p-u*.15,3.1)
        lid=lid+rod(p+u*.15,p+u*3.5,3.1)
        hole=rod(p-u*10,p+u*10,1.1);rear=rear-hole;lid=lid-hole
        head=p+u*3.8
        hardware[label+'_M2x10']=stepped_bolt(head,u,v,w) if DISTAL_LOAD else rod(head,head-u*10,1.)+rod(head,head+u*1.6,1.9)
        hardware[label+'_head_washer']=rod(p+u*3.5,head,2.5)-hole
        hardware[label+'_back_washer']=rod(p-u*3.5,p-u*3.8,2.5)-hole
        hardware[label+'_nut']=mf.Manifold.cylinder(1.6,4/np.sqrt(3),4/np.sqrt(3),6).transform(np.c_[np.stack([v,-w,-u],1),p-u*3.8])-rod(p-u*3,p-u*7,1.05)
    rear=(rear-cavity-window).as_original().simplify(.01)
    # Original cap may extend into the split-front region; relieve only the
    # retainer's actual envelope, preserving source M3 cap/pad attachment.
    clearance=mf.Manifold.batch_boolean([lid]+[lid.translate(d*s*.12) for d in np.eye(3) for s in [-1,1]],mf.OpType.Add)
    rear=(rear-clearance).as_original().simplify(.01);lid=lid.as_original().simplify(.01);unstopped=rear
    if DISTAL_LOAD:
        # Positive printed shoulders carry bolt clamp load and preserve the
        # 0.3mm shell split. Without them tightening could shift the lid and
        # consume the ball clearance; a nominal open gap is not a hard stop.
        for p in eyes:rear=rear+(rod(p-u*.15,p+u*.15,2.5)-rod(p-u,p+u,1.1))
        rear=rear.as_original().simplify(.01)
    terminal=(mf.Manifold.sphere(3.5,32).translate(C)+rod(C+u*3,C+u*6.2,1.2)+rod(C+u*6,C+u*9,1.6))-rod(C+u*3.3,C+u*10,.45)
    parts={'distal_socket_base':rear,'socket_retainer':lid,'ball_crimp_blank':terminal.as_original().simplify(.005),**hardware}
    return parts,eyes,window,unstopped

def run():
    DEST.mkdir(exist_ok=True)
    assembly=json.loads((BASE/'assembly.json').read_text(encoding='utf8'));lever=json.loads((BASE/'design.json').read_text(encoding='utf8'));rows={r['key']:r for r in assembly['instances']}
    baseline={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()}
    placements={p['id']:np.array(p['translation_display_mm']) for p in assembly['placements']};all_reports=[]
    for fid,key in [('FJ1486','flexor'),('FJ1477','extensor')]:
        folder=DEST/fid;folder.mkdir(exist_ok=True);design=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));F=np.array(design['basis_columns_display']);pad=np.array(design['pad_centres_display_mm'][1])+placements[fid];capkey=fid+'_distal_rigid_cap';cap=baseline[capkey]
        target=np.array(lever['eyelets_display_mm'][key]);side=1 if (target-pad)@F[:,0]>=0 else -1
        candidates=[];selected=None;best=None
        for position in [(14.5,8.,8.),(16.,9.,9.),(17.5,10.,10.)]:
            C=pad+F@np.array([side*position[0],position[1],position[2]]);nodes=np.array(design['wire_polyline_display_mm']);active_pull=unit(nodes[-2]-nodes[-1]);u=-active_pull if DISTAL_LOAD else unit(target-C);zero=F[:,2]-u*(F[:,2]@u)
            if np.linalg.norm(zero)<.1:zero=F[:,0]-u*(F[:,0]@u)
            zero=unit(zero);other=np.cross(u,zero);root_local=np.array([side*7,7,5. if DISTAL_LOAD else 4.3]);root=pad+F@root_local
            for roll in [0,45,90,135]:
                v=zero*np.cos(np.radians(roll))+other*np.sin(np.radians(roll));parts,eyes,window,unstopped=construct(cap,C,u,v,root)
                connected={name:len([s for s in p.decompose() if s.volume()>.001]) for name,p in parts.items()};hard=[];soft=[]
                for (a,sa),(b,sb) in itertools.combinations(parts.items(),2):
                    volume=intersection(sa,sb)
                    if volume>.001:hard.append({'pair':[a,b],'volume_mm3':volume})
                for a,sa in parts.items():
                    for b,sb in baseline.items():
                        if b==capkey:continue
                        volume=intersection(sa,sb)
                        if volume>.001:(soft if rows[b]['kind']=='TPU' else hard).append({'pair':[a,b],'volume_mm3':volume})
                row={'position_local_mm':[side*position[0],position[1],position[2]],'bridge_root_local_mm':root_local.tolist(),'retainer_roll_deg':roll,'connected_components':connected,'hard_interferences':hard,'TPU_contact_review_flags':soft,'ball_centre_display_mm':C.tolist(),'neutral_terminal_axis':u.tolist(),'retainer_v_axis':v.tolist(),'M2_centres_display_mm':[p.tolist() for p in eyes]}
                candidates.append(row);score=(sum(v!=1 for v in connected.values()),len(hard),sum(h['volume_mm3'] for h in hard),len(soft),np.linalg.norm(np.array(position)-[11.3,4,4.15]))
                if best is None or score<best:best=score;selected=(row,parts,unstopped)
                print(fid,position,roll,'components',set(connected.values()),'hard',len(hard),'soft',len(soft),flush=True)
            if best[0]==best[1]==0 and best[3]==0:break
        row,parts,unstopped=selected;exports=[write(folder/(name+'.stl'),solid) for name,solid in parts.items()]
        if DISTAL_LOAD:
            source_folder=folder/'native_base_source';source_folder.mkdir(exist_ok=True);source_part=write(source_folder/'distal_socket_base.stl',unstopped)
            (source_folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[source_part]}),encoding='utf8')
        report={'stage':'M5 physical distal spherical cable terminal candidate','muscle':fid,'selected':row,'candidate_comparison':candidates,'parts':{r['name']:r['components'] for r in exports},'ball_radius_mm':3.5,'socket_inner_radius_mm':3.7,'socket_outer_radius_mm':5.2,'nominal_articulation_cone_half_angle_deg':25,'neck_radius_mm':1.2,'crimp_body_radius_mm':1.6,'uncrimped_wire_bore_radius_mm':.45,'nominal_wire_radius_mm':.3,'split_half_gap_mm':.3,'fasteners':'two nominal M2x10 bolts, two nuts, four washers; threads/preload absent','source_cap_key':capkey,'source_cap_OBJ_sha256':hashlib.sha256((ROOT/rows[capkey]['file']).read_bytes()).hexdigest(),'source_assembly_sha256':hashlib.sha256((BASE/'assembly.json').read_bytes()).hexdigest(),'force_path':'TPU distal pad -> existing M3 end-cap attachment -> integral socket bridge -> spherical ball head -> crimped passive tendon -> later real guide/forearm attachment. Crimp/strength/guide route not verified.','pending':['retainer/ball/fastener installation','articulation vs actual components and retention','CAD roundtrip and printability','passive steel termination/grip, supported cable route','TPU deformation and full-arm installed motion','strength, tolerance, wear and loads'],'frozen':False}
        report.update(neutral_outlet_mode='distal_load_opposes_active_pull' if DISTAL_LOAD else 'direct_to_old_lever_reference',neutral_active_wire_pull_direction=active_pull.tolist(),configured_neutral_outlet_opposition=-float(np.array(row['neutral_terminal_axis'])@active_pull))
        report['integral_clamp_stops']={'present':DISTAL_LOAD,'count':2 if DISTAL_LOAD else 0,'annular_outer_radius_mm':2.5,'bore_radius_mm':1.1,'height_mm':.3,'purpose':'rear-ear shoulders contact retainer faces to set shell split under nominal tightening; load/creep not validated'}
        (folder/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');(folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':exports}),encoding='utf8');all_reports.append({'muscle':fid,'selected':row})
    (DEST/'comparison.json').write_text(json.dumps({'scope':'real split socket endpoints, not completed tendon routing or TPU drive','modules':all_reports},indent=2),encoding='utf8')
if __name__=='__main__':run()
