"""Real M10-R1 cap ports, steel envelopes and extracted manufactured guide pieces.

Interior guide centres are free numerical variables. Their moved poses are a
geometric hypothesis, not a solved TPU body or imaginary fixed supports. Report
each scope separately: bones, metal, short guide lumen, and unknown full lattice.
"""
from build_direct_muscle_yoke import *
from study_extensor_bone_corridor import BoneDistance
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from pathlib import Path

def capsule_chain(points,radius):
    return mf.Manifold.batch_boolean([rod(p,q,radius,24) for p,q in zip(points,points[1:])]+[mf.Manifold.sphere(radius,16).translate(p) for p in points[1:-1]],mf.OpType.Add)

class NearbyDistance:
    def __init__(self,solid):
        self.distance=BoneDistance(solid);self.box=np.array(solid.bounding_box())
        if len(self.distance.tri)<1024:
            class AllTriangles:
                def __init__(self,count):self.count=count
                def query(self,points,k):return None,np.broadcast_to(np.arange(self.count),(len(points),self.count))
            self.distance.tree=AllTriangles(len(self.distance.tri))
    def query(self,points):
        delta=points-np.clip(points,self.box[:3],self.box[3:]);value=np.linalg.norm(delta,axis=1);gradient=delta/np.maximum(value[:,None],1e-12);near=value<5
        if near.any():value[near],gradient[near]=self.distance.query(points[near])
        return value,gradient

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';candidate=parent/'captive_keepers';dest=candidate/'wire_pose_study';dest.mkdir(exist_ok=True);a=json.loads((candidate/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if r['kind']!='TPU'};motion=json.loads((parent/'conditional_motion.json').read_text())['poses'];j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);bone_keys=['FJ3368','FJ3391'];distance={k:BoneDistance(world[k]) for k in bone_keys};inputs=[candidate/'assembly.json',parent/'conditional_motion.json']+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};results=[];guide_sources=[]
    metal_aware='--metal-aware' in sys.argv
    if metal_aware:dest=candidate/'wire_pose_metal_aware';dest.mkdir(exist_ok=True)
    regularized='--regularized' in sys.argv;normal_bow='--normal-bow' in sys.argv;critical='--critical-flexor' in sys.argv
    if regularized or normal_bow:dest=candidate/('wire_pose_regularized_critical' if critical else 'wire_pose_regularized');dest.mkdir(exist_ok=True)
    if '--fine-wire' in sys.argv:dest=candidate/'wire_pose_regularized_fine';dest.mkdir(exist_ok=True)
    separated='--separate-guides' in sys.argv;continuation='--continuation' in sys.argv
    if separated:dest=candidate/'wire_pose_separated_guides';dest.mkdir(exist_ok=True)
    metal_queries={k:NearbyDistance(s) for k,s in world.items() if metal_aware and (k=='paired_bone_frame' or '_keeper' in k)}
    angles=[0,30,60,90,120,132,140,150] if '--quick' in sys.argv else sorted(set(p['elbow_deg'] for p in motion))
    if critical:angles=[120,132,140,150]
    if '--single-critical' in sys.argv:angles=[132]
    for fid in (['FJ1486'] if critical else ['FJ1486','FJ1477']):
        folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());hashes[str(folder/'design.json')]=hashlib.sha256((folder/'design.json').read_bytes()).hexdigest();offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));A,B=np.array(d['pad_centres_display_mm'])+offset;z=np.array(d['basis_columns_display'])[:,2];guides=d['actual_short_tube_guides'];centres=np.array([g['centre_display_mm'] for g in guides])+offset;rest=np.array([A-z*5.15,A-z*2.15,*centres,B+z*2.15,B+z*4.6]);N=len(rest);s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(rest,axis=0),axis=1))];s/=s[-1];weights=list(np.eye(N))
        for i in range(N-1):
            for t in (np.arange(1,32)/32 if '--fine-wire' in sys.argv else np.arange(1,8)/8):
                w=np.zeros(N);w[i]=1-t;w[i+1]=t;weights.append(w)
        weights=np.array(weights);tube_source=[];muscle=load(folder/'TPU_muscle.obj')[0].translate(offset)
        feature=next(f for f in json.loads((candidate/'design.json').read_text())['features'] if f['muscle']==fid);screen_path=DEST/('cap_swivel_flexor' if fid=='FJ1486' else 'cap_swivel_study')/'screen.json';screen=json.loads(screen_path.read_text());e=np.array(screen['trunnion_axis_zero_display']);bow_normal=unit(np.cross(np.array(feature['keeper_bolt_axis_display']),e));hashes[str(candidate/'design.json')]=hashlib.sha256((candidate/'design.json').read_bytes()).hexdigest();hashes[str(screen_path)]=hashlib.sha256(screen_path.read_bytes()).hexdigest()
        for g,c in zip(guides,centres):
            tangent=np.array(g['axis_display']);half=g['length_mm']/2;outer=g['outer_radius_mm'];tube=muscle^rod(c-tangent*half,c+tangent*half,outer+.015,48);tube_source.append(tube)
        nominal= capsule_chain(rest,.3);nominal_hits=[float(intersection(nominal,t)) for t in tube_source];guide_sources.append({'muscle':fid,'tube_count':len(tube_source),'extraction':'actual manufactured TPU intersect actual nominal guide-cylinder region, not replacement cylinders','source_neutral_wire_max_guide_overlap_mm3':max(nominal_hits),'extract_volumes_mm3':[float(t.volume()) for t in tube_source]});print(fid,'extracted real guide pieces',len(tube_source),'baseline maxsteelhit',max(nominal_hits),flush=True);last=None;last_forearm=None
        pair_i,pair_j=np.array(list(itertools.combinations(range(len(guides)),2))).T;guide_bounds=np.array([np.hypot(g['outer_radius_mm']+.015,g['length_mm']/2) for g in guides]);required_distance=guide_bounds[pair_i]+guide_bounds[pair_j]+.1
        for angle in angles:
            pose=next(p for p in motion if p['elbow_deg']==angle);Tf=np.array(pose['forearm_pose']);Tc=np.array(pose['cap_poses'][fid]);ends=rest[-2:]@Tc[:,:3].T+Tc[:,3];initial=rest+(ends[-1]-rest[-1])*s[:,None]
            if last is not None:initial=last+(ends[-1]-last[-1])*s[:,None]
            if continuation and last is not None:
                delta_R=Tf[:,:3]@last_forearm[:,:3].T;initial=last+s[:,None]*((last-E)@delta_R.T+E-last)
            initial[:2]=rest[:2];initial[-2:]=ends
            if normal_bow and (not continuation or last is None) and (not continuation or angle>=110):initial[2:-2]+=np.sin(np.pi*s[2:-2])[:,None]*(Tf[:,:3]@bow_normal)*12
            u=s[2:-2];w=6*u**5-15*u**4+10*u**3;rest_target=centres+w[:,None]*(centres@Tc[:,:3].T+Tc[:,3]-centres)
            def unpack(x):return np.vstack([rest[:2],x.reshape(-1,3),ends])
            def objective(x):
                P=unpack(x);links=np.diff(P,axis=0);length=np.sqrt((links*links).sum(1)+1e-12);curv=np.diff(P,n=2,axis=0);g=np.zeros_like(P);u=links/length[:,None];g[:-1]-=u;g[1:]+=u;g[:-2]+=.03*curv;g[1:-1]-=.06*curv;g[2:]+=.03*curv;cost=float(length.sum()+.015*(curv*curv).sum())
                if regularized:
                    deviation=P[2:-2]-rest_target;cost+=float(.05*(deviation*deviation).sum());g[2:-2]+=.1*deviation
                return cost,g[2:-2].reshape(-1)
            def constraints(x):
                P=unpack(x);Q=weights@P;values=[];jac=[]
                for key,query in distance.items():
                    local=(Q-Tf[:,3])@Tf[:,:3] if rows[key]['owner']=='elbow' else Q;value,g=query.query(local)
                    if rows[key]['owner']=='elbow':g=g@Tf[:,:3].T
                    J=weights[:,:,None]*g[:,None,:];values.extend(value-2.7-.3);jac.append(J[:,2:-2].reshape(len(Q),-1))
                for key,query in metal_queries.items():
                    local=(Q-Tf[:,3])@Tf[:,:3];value,g=query.query(local);g=g@Tf[:,:3].T;J=weights[:,:,None]*g[:,None,:];values.extend(value-.45);jac.append(J[:,2:-2].reshape(len(Q),-1))
                    local_centres=(P[2:-2]-Tf[:,3])@Tf[:,:3];value,g=query.query(local_centres);g=g@Tf[:,:3].T;J=np.zeros((len(value),len(value),3));J[np.arange(len(value)),np.arange(len(value))]=g;values.extend(value-2.8);jac.append(J.reshape(len(value),-1))
                if separated:
                    C=P[2:-2];delta=C[pair_i]-C[pair_j];length=np.maximum(np.linalg.norm(delta,axis=1),1e-12);g=delta/length[:,None];J=np.zeros((len(pair_i),len(C),3));J[np.arange(len(pair_i)),pair_i]=g;J[np.arange(len(pair_i)),pair_j]=-g;values.extend(length-required_distance);jac.append(J.reshape(len(pair_i),-1))
                return np.array(values),np.vstack(jac)
            solution=minimize(lambda x:objective(x)[0],initial[2:-2].reshape(-1),jac=lambda x:objective(x)[1],constraints=[{'type':'ineq','fun':lambda x:constraints(x)[0],'jac':lambda x:constraints(x)[1]}],method='SLSQP',options={'maxiter':180,'ftol':1e-7});P=unpack(solution.x);wire=capsule_chain(P,.3);corridor=capsule_chain(P,2.7);posed={}
            for key,body in world.items():
                owner=rows[key]['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;posed[key]=body.transform(T) if T is not None else body
            steel_hits=[];bone_hits=[]
            for key,body in posed.items():
                v=intersection(wire,body)
                if v>.001:steel_hits.append({'part':key,'volume_mm3':v,'nominal_uncrimped_stop':rows[key].get('nominal_packing_only',False)})
                if key in bone_keys:
                    v=intersection(corridor,body)
                    if v>.001:bone_hits.append({'bone':key,'volume_mm3':v})
            tube_checks=[]
            for i,(guide,c,tube) in enumerate(zip(guides,centres,tube_source)):
                index=i+2;incoming=unit(P[index]-P[index-1]);outgoing=unit(P[index+1]-P[index]);tangent=unit(incoming+outgoing);R=Rotation.align_vectors([tangent],[np.array(guide['axis_display'])])[0].as_matrix();T=np.c_[R,P[index]-R@c];moved=tube.transform(T);v=intersection(wire,moved);collisions=[]
                for key,body in posed.items():
                    overlap=intersection(moved,body)
                    if overlap>.001:collisions.append({'part':key,'volume_mm3':overlap})
                tube_checks.append({'guide':i,'pose_affine_mm':T.tolist(),'incoming_outgoing_turn_deg':float(np.degrees(np.arccos(np.clip(incoming@outgoing,-1,1)))),'steel_material_intersection_mm3':v,'rigid_intersections':collisions,'actual_TPU_embedded_guide_pose_verified':False})
            separation_margin=float((np.linalg.norm(P[2:-2][pair_i]-P[2:-2][pair_j],axis=1)-required_distance).min());accepted=bool(solution.success and not bone_hits and not steel_hits and not any(t['steel_material_intersection_mm3']>.001 or t['rigid_intersections'] for t in tube_checks) and (not separated or separation_margin>=-1e-6))
            record={'muscle':fid,'elbow_deg':angle,'optimizer_success':bool(solution.success),'optimizer_message':str(solution.message),'optimizer_minimum_constraint_mm':float(constraints(solution.x)[0].min()),'centreline_display_mm':P.tolist(),'geometric_wire_length_mm':float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum()),'bone_corridor_intersections':bone_hits,'steel_rigid_intersections':steel_hits,'guide_checks':tube_checks,'minimum_enclosing_sphere_separation_margin_mm':separation_margin if separated else None,'maximum_guide_centre_shift_from_previous_accepted_sample_mm':float(np.linalg.norm(P[2:-2]-last[2:-2],axis=1).max()) if last is not None else None,'local_proxy_geometry_acceptance':accepted,'actual_TPU_body_deformation_verified':False,'actual_friction_or_drive_verified':False};results.append(record);last=P if accepted else None;last_forearm=Tf if accepted else None
            print(fid,'q',angle,'opt',solution.success,'bone',len(bone_hits),'steelhard',len(steel_hits),'guidewire',sum(t['steel_material_intersection_mm3']>.001 for t in tube_checks),'guiderigid',sum(bool(t['rigid_intersections']) for t in tube_checks),flush=True)
            (dest/'partial_checkpoint.json').write_text(json.dumps({'scope':'in-progress geometric hypothesis only, not actual TPU or drive proof','poses':results,'guides':guide_sources},indent=2),encoding='utf8')
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());report={'scope':'actual swivel-cap port endpoints, numerical free-centreline path,0.6mmsteel and rigidly posed source-manufactured short-guide pieces. These are hypothetical TPU-guide poses, NOT lattice deformation, equilibrium, friction ordrive proof; no actual model guides are fixed to virtual supports','quick':'--quick' in sys.argv,'critical_flexor_only':critical,'single_critical_case':'--single-critical' in sys.argv,'sampling_per_link':32 if '--fine-wire' in sys.argv else 8,'guide_separation_constraint':'enclosing spheres for actual guide-cylinder clipping region plus0.1mm; conservative rigid-piece compatibility assumption, not calibrated TPU stiffness' if separated else None,'continuation_predictor':continuation,'small_metal_distance_method':'all triangles for solids<1024faces; nearest128 forlarger solids; independent actualsolid Booleans finalgate','regularizer':'0.05mm^-1 squared-distance to smooth rest/cap interpolation; numerical shape regularizer, not calibrated TPU elasticity or external anchors' if regularized else None,'normal_bow_initial_guess_mm':12 if normal_bow else 0,'metal_aware':metal_aware,'constrained_metal_parts':list(metal_queries),'source_sha256':hashes,'guide_sources':guide_sources,'poses':results,'frozen':False};(dest/('wire_guides_quick.json' if '--quick' in sys.argv else 'wire_guides_all.json')).write_text(json.dumps(report,indent=2),encoding='utf8')

if __name__=='__main__':run()
