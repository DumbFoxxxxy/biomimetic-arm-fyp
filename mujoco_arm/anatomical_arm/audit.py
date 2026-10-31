"""Exact mesh bone audit; hardware and TPU receive independent neutral checks.
Kinematic sweeps are not actuation or compliance predictions.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build import *
import mujoco

def load(n):
    v,f=read_stl(str(OUT/'assets'/(n+'.stl')))
    return mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))

def run():
    m=_load_project_model(str(OUT/'arm.xml'));d=mujoco.MjData(m)
    frames=json.loads((OUT/'body_frames.json').read_text(encoding='utf8'));orig={k:np.array(v) for k,v in frames['origins_mm'].items()};owners=frames['bone_owners'];base={n:load(n) for n in owners}
    bids={n:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,n) for n in orig}
    jids={mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,i):int(m.jnt_qposadr[i]) for i in range(m.njnt)}
    muscles=json.loads((OUT/'muscle_parameters.json').read_text(encoding='utf8'));rows=[];posecache={};paths={}
    pairs=[(a,b) for a,b in itertools.combinations(base,2) if owners[a]!=owners[b]]
    def pose(q):
        key=tuple(np.round(q,9))
        if key in posecache:return posecache[key]
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);world={};transforms={}
        for n,s in base.items():
            owner=owners[n];i=bids[owner];R=d.xmat[i].reshape(3,3);t=d.xpos[i]*1000-R@orig[owner];T=np.c_[R,t];world[n]=s.transform(T);transforms[owner]=T
        hits=[];gap=5.
        for a,b in pairs:
            aa=np.array(world[a].bounding_box());bb=np.array(world[b].bounding_box());delta=np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:]);lower=np.linalg.norm(np.maximum(delta,0))
            if lower<5:
                volume=float((world[a]^world[b]).volume()) if np.all(delta<0) else 0.
                if volume>1e-3:hits.append({'pair':[a,b],'volume_mm3':volume})
                gap=min(gap,float(world[a].min_gap(world[b],5)))
        E=np.array(json.loads((OUT/'source_audit.json').read_text(encoding='utf8'))['landmarks_mm']['E']);GH=orig['humerus'];T=transforms['humerus'];direction=T[:,:3]@(E-GH);zero=E-GH
        elev=float(np.degrees(np.arccos(np.clip(unit(direction)@unit(zero),-1,1))))
        R=d.xmat[bids['scapula']].reshape(3,3);v=R@np.array([1.,0,0]);scap=float(np.degrees(np.arctan2(v[2],v[0])))
        tendon=[]
        for md in muscles:
            C=np.array(md['C_wire_terminal_and_tendon_mm']);B=np.array(md['B_bone_insertion_mm']);Tc=transforms['humerus'];Tb=transforms['forearm'];C=Tc[:,:3]@C+Tc[:,3];B=Tb[:,:3]@B+Tb[:,3]
            tendon.append(float(np.linalg.norm(B-C)))
        value={'hits':hits,'minimum_bone_gap_mm_capped_5':gap,'humerus_elevation_deg':elev,'scapular_frontal_rotation_deg':scap,'neutral_C_to_moving_B_distance_mm':tendon}
        posecache[key]=value;return value
    actions=[('K1',{'elbow':145},145,116),('K2',{'gh_flex':120},120,96),('K3',{'gh_abd':120},120,96),('K4',{'gh_rot':70},70,56),('K5',{'gh_rot':-90},90,72),('K6',{'ac_up':50},50,36),('K7',{'ac_up':50,'gh_abd':120},170,132),('SCAP_PROTRACT',{'sc_protract':20},20,None),('SCAP_RETRACT',{'sc_protract':-20},20,None),('SCAP_COMBINED_UP',{'sc_elev':20,'ac_up':30},50,None),('V2_SLIDE_FWD_COORDINATE',{'sc_protract':5},5,None),('V2_SLIDE_BACK_COORDINATE',{'sc_protract':-5},5,None),('PECTORAL',{'gh_flex':25,'gh_abd':8,'gh_rot':-20},25,None)]
    for aid,target,goal,floor in actions:
        end=np.zeros(m.nq)
        for n,val in target.items():end[jids[n]]=val
        steps=int(np.ceil(max(abs(end))));forward=[];first=None
        for i in range(steps+1):
            f=i/steps;q=end*f;p=pose(q);forward.append({'fraction':f,'q_deg':q.tolist(),**p})
            if p['hits'] and first is None:first=i
        bound=None;actual=0
        if first is None:actual=1.
        elif first>0:
            lo=(first-1)/steps;hi=first/steps
            while max(abs(end))*(hi-lo)>.1:
                mid=(lo+hi)/2
                if pose(end*mid)['hits']:hi=mid
                else:lo=mid
            actual=lo;bound=pose(end*hi)
        else:bound=forward[0]
        # Explicitly replay return samples rather than assuming valid return.
        return_hits=0
        for p in reversed(forward):
            q=np.array(p['q_deg']);posecache.pop(tuple(np.round(q,9)),None);rp=pose(q);return_hits+=bool(rp['hits'])
        accepted=pose(end*actual);actualmetric=accepted['humerus_elevation_deg'] if aid=='K7' else abs(goal)*actual
        row={'action':aid,'target':goal,'floor':floor,'commanded_joint_endpoint_deg':target,'actual_continuous_sampled_prefix':actualmetric,'metric':'actual humerus direction from zero' if aid=='K7' else 'prescribed coordinate/path parameter, not total anatomical ROM','first_limit':bound['hits'] if bound else [],'minimum_bone_gap_mm_capped_5':min(p['minimum_bone_gap_mm_capped_5'] for p in forward[:(first if first else len(forward))]),'actual_humerus_angle_deg':accepted['humerus_elevation_deg'],'scapular_frontal_rotation_deg':accepted['scapular_frontal_rotation_deg'],'forward_samples':len(forward),'return_samples':len(forward),'return_colliding_samples':return_hits,'bone_floor_pass':None if floor is None else actualmetric>=floor,'hardware_pass':False,'actuation_pass':None}
        rows.append(row);paths[aid]=forward;print(aid,round(actualmetric,2),row['first_limit'],flush=True)
    neutral=pose(np.zeros(m.nq));hardware=json.loads((OUT/'hardware.json').read_text(encoding='utf8'));hw=[]
    for part in hardware['parts']:
        n=part['id'];s=load(n)
        for bone,t in base.items():
            vol=(s^t).volume()
            if vol>.001:hw.append({'hardware':n,'bone':bone,'overlap_mm3':vol,'classification':'unresolved; pocket/support connection not designed'})
    muscle_checks=[]
    for md in muscles:
        s=load(md['id']);hits=[]
        for bone,t in base.items():
            vol=(s^t).volume()
            if vol>.001:hits.append({'bone':bone,'volume_mm3':vol,'classification':'neutral muscle penetration requires rerouting; not accepted compliant contact'})
        muscle_checks.append({'id':md['id'],'closed':str(s.status()),'material_components':len([x for x in s.decompose() if x.volume()>0]),'bone_overlaps':hits})
    dump(OUT/'range_audit.json',{'status':'bone-only audit; candidate not frozen','overlap_volume_reporting_threshold_mm3':.001,'angular_step_max_deg':1,'boundary_resolution_deg':.1,'gap_search_cap_mm':5,'all_adjacent_bone_pairs_checked':True,'legal_contacts':[],'radius_ulna_neutral_gap_mm_capped_5':base['radius'].min_gap(base['ulna'],5),'same_rigid_body_pair':'radius-ulna neutral gap checked separately; no relative motion','unsupported_actions':['independent scapular internal/external rotation absent in AC1 chain','translational scapular slide not an independent coordinate; needs linkage IK'],'neutral':neutral,'actions':rows,'hardware_neutral_overlaps':hw,'muscle_neutral':muscle_checks,'acceptance':'No assembly or drive PASS while hardware attachment/clearance unresolved.'})
    dump(OUT/'motion_paths.json',paths)
    with (OUT/'acceptance.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['action','target','floor','actual_continuous_sampled_prefix','bone_floor_pass','hardware_pass','first_limit']);writer.writeheader();writer.writerows({k:r[k] for k in writer.fieldnames} for r in rows)
    print('hardware overlaps',len(hw),'muscle checks',muscle_checks,flush=True)

if __name__=='__main__':run()
