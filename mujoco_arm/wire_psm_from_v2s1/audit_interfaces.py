"""Dense actual swivel paths, belly-space scenario and finite sheave envelopes."""
from common import *
import manifold3d as md
import itertools

def segment_gap(a,b,c,d):
    u=b-a;v=d-c;r=a-c;aa=u@u;ee=v@v;f=v@r
    if aa<1e-16 and ee<1e-16:return float(np.linalg.norm(a-c))
    if aa<1e-16:s=0.;t=np.clip(f/ee,0,1)
    else:
        cc=u@r
        if ee<1e-16:t=0.;s=np.clip(-cc/aa,0,1)
        else:
            bb=u@v;den=aa*ee-bb*bb;s=np.clip((bb*f-cc*ee)/den,0,1) if den>1e-16 else 0.;t=(bb*s+f)/ee
            if t<0:t=0.;s=np.clip(-cc/aa,0,1)
            elif t>1:t=1.;s=np.clip((bb-cc)/aa,0,1)
    return float(np.linalg.norm(a+s*u-c-t*v))

def wire_pair_gap(a,b,threshold=2.):
    if np.any(a.min(0)-threshold>b.max(0)) or np.any(b.min(0)-threshold>a.max(0)):return 1e6
    minimum=1e6
    for x,y in zip(a[:-1],a[1:]):
        mask=np.all(np.minimum(b[:-1],b[1:])<=np.maximum(x,y)+threshold,axis=1)&np.all(np.maximum(b[:-1],b[1:])>=np.minimum(x,y)-threshold,axis=1)
        for c,d in zip(b[:-1][mask],b[1:][mask]):minimum=min(minimum,segment_gap(x,y,c,d)-1.)
    return minimum

def run():
    spec=json.loads((OUT/'swivel_interfaces.json').read_text());modules=spec['modules'];z=np.load(OUT/'swivel_interfaces.npz');m,d=compile_tree(ET.parse(OUT/'psm_module_candidate.xml'));parts=solids(m,d)
    cb=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,x['carrier']) for x in modules];tb=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,x['insertion_body']) for x in modules];A=z['A'];G=z['G'];B=z['B'];axis=z['axis'];cache={};worst={};actions={};rows=[];zero=[]
    normals=[]
    for i,x in enumerate(modules):
        n=np.array(x['belly_normal_local_mj']);n-=axis[i]*(n@axis[i]);n/=max(np.linalg.norm(n),1e-12);v=np.cross(n,axis[i]);normals.append((v,n))
    def overlap(shape,part,key):
        if key in cache:return cache[key]
        box=np.array(shape.bounding_box())
        if np.any(box[:3]>=part['hi']) or np.any(box[3:]<=part['lo']):vol=0.
        else:vol=float((shape^part['solid']).volume())
        cache[key]=vol;return vol
    for action,(start,end) in spec['action_slices'].items():
        first_hard=None;first_belly=None;pair_first=None
        for k in range(start,end):
            q=z['q_deg'][k];d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);R=d.xmat.reshape(-1,3,3);P=d.xpos*1000;wires=[];issues=[];belly_issues=[]
            for i,module in enumerate(modules):
                carrier=cb[i];a=A[i]@R[carrier].T+P[carrier];g=G[i]@R[carrier].T+P[carrier];b=B[i]@R[tb[i]].T+P[tb[i]];r=swivel_route(a,g,b)
                if r is None or r['wrap_deg']>180+1e-6:issues.append({'module':module['id'],'kind':'guide_geometry','wrap_deg':r['wrap_deg'] if r else None});wires.append(None);continue
                wires.append(r['path'])
                if not z['feasible'][k,i]:issues.append({'module':module['id'],'kind':'passive_tendon_C_travel'})
                s=float(np.clip(z['s_mm'][k,i],20,module['available_AG_mm']-24));v,n=normals[i];belly=md.Manifold.cube([s,34,14],center=True).transform(np.column_stack([axis[i],v,n,A[i]+axis[i]*s/2]))
                wheel=md.Manifold.cylinder(6,9.5,circular_segments=24,center=True).transform(np.column_stack([(g-a)/np.linalg.norm(g-a),np.cross(r['axis'],(g-a)/np.linalg.norm(g-a)),r['axis'],r['center']]))
                for name,part in parts.items():
                    body=part['body'];local=(r['path']-P[body])@R[body];key=('wire',i,name,np.round(local,4).tobytes());gap=cache.get(key)
                    if gap is None:gap=local_gap(local,part);cache[key]=gap
                    if gap<2:issues.append({'module':module['id'],'kind':'wire_bone','bone':name,'gap_mm':float(gap)})
                    # Complete outer envelopes, no joint-pair exceptions.
                    rel=R[body].T@R[carrier];shift=R[body].T@(P[carrier]-P[body]);localbelly=belly.transform(np.column_stack([rel,shift]));key=('belly',i,name,round(s,3),np.round(np.r_[rel.ravel(),shift],4).tobytes());vol=overlap(localbelly,part,key)
                    if vol>.001:belly_issues.append({'module':module['id'],'kind':'TPU_14x34_envelope_risk','bone':name,'volume_mm3':vol})
                    localwheel=wheel.transform(np.column_stack([R[body].T,-R[body].T@P[body]]));key=('guide',i,name,np.round(np.r_[R[body].T@(r['center']-P[body]),R[body].T@r['axis'],R[body].T@((g-a)/np.linalg.norm(g-a))],4).tobytes());vol=overlap(localwheel,part,key)
                    if vol>.001:issues.append({'module':module['id'],'kind':'finite_guide_bone','bone':name,'volume_mm3':vol})
            pairs=[]
            for i,j in itertools.combinations(range(16),2):
                if wires[i] is None or wires[j] is None:continue
                gap=wire_pair_gap(wires[i],wires[j])
                if gap<1:pairs.append({'modules':[modules[i]['id'],modules[j]['id']],'surface_gap_mm':float(gap)})
            if first_hard is None and issues:first_hard={'index':k-start,'q_deg':q.tolist(),'issues':issues}
            if first_belly is None and belly_issues:first_belly={'index':k-start,'q_deg':q.tolist(),'issues':belly_issues}
            if pair_first is None and pairs:pair_first={'index':k-start,'q_deg':q.tolist(),'pairs':pairs}
            for issue in issues+belly_issues:
                key=(issue['module'],issue['kind'],issue.get('bone'));worst.setdefault(key,{**issue,'action':action,'index':k-start})
            if k==0:zero={'hard_issues':issues,'TPU_envelope_risks':belly_issues,'wire_pair_issues':pairs}
        row={'action':action,'sample_count':end-start,'first_hard_geometry_failure':first_hard,'first_TPU_envelope_risk':first_belly,'first_wire_pair_failure':pair_first};rows.append(row);print('interfaces',action,'first hard',None if first_hard is None else first_hard['index'],'belly',None if first_belly is None else first_belly['index'],'pairs',None if pair_first is None else pair_first['index'],flush=True)
        report={'scope':'807 frozen bone poses; actual R8 swivelling sheave arcs; 1mm steel with2mm bone clearance and1mm inter-wire surface clearance; finite guide outer envelope radius9.5 width6; TPU14x34 rectangular occupied-space scenario, not calibrated deformed lattices; bracket roots/hardware/motor feeds not certified','zero':zero,'actions':rows,'failure_groups':list(worst.values()),'frozen':False,'hardware_ready':False,'TPU_body_packing_verified':False}
        (OUT/'interfaces_geometry_audit.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':run()
