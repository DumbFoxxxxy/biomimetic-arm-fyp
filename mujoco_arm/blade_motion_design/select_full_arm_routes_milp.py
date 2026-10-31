"""Jointly select collision-filtered fairleads with pull-only closure at each pose."""
import json,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import lil_matrix,coo_matrix
from shoulder_model import EXPORT,inputs,SOURCE
import sys
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances
from arm_targets import target_actions
from shoulder_model import v8

def run():
    cfg,fit,existing=inputs();pools=json.loads((EXPORT/'full_arm_route_pool.json').read_text());arrays=np.load(EXPORT/'full_arm_route_pool.npz');weights=np.concatenate([arrays[f'g{g}_w'] for g in range(8)]);sizes=[len(p) for p in pools];starts=np.r_[0,np.cumsum(sizes)];N=len(weights);indices=[0,5,11,17,23,29,35,77,83];P=len(indices);V=N*(P+1)
    # X picks the fixed fairlead. Y gives the nonnegative normalized tension
    # of that selected line independently at each target pose.
    rows=[];cols=[];vals=[];lower=[];upper=[]
    def add(coeff,lo,hi):
        r=len(lower)
        for c,v in coeff:rows.append(r);cols.append(int(c));vals.append(float(v))
        lower.append(lo);upper.append(hi)
    for g in range(8):add([(j,1) for j in range(starts[g],starts[g+1])],1 if g<4 else 0,2)
    add([(j,1) for j in range(N)],8,8)
    for k,pose in enumerate(indices):
        base=N*(k+1)
        for axis in range(5):add([(base+j,weights[j,pose,axis]) for j in range(N)],0,0)
        add([(base+j,1) for j in range(N)],1,1)
        for g in range(4):add([(base+j,1) for j in range(starts[g],starts[g+1])],1e-4,np.inf)
        for j in range(N):
            add([(base+j,1),(j,-1)],-np.inf,0)
            add([(base+j,1),(j,-1e-5)],0,np.inf)
    objective=np.zeros(V)
    intg=np.r_[np.ones(N),np.zeros(N*P)];bounds=Bounds(np.zeros(V),np.ones(V));cuts=set()
    poses=[];seen=set()
    for name,seq in target_actions().items():
        for a,b in zip(seq[:-1],seq[1:]):
            for f in np.linspace(0,1,5,endpoint=False):
                z=a*(1-f)+b*f;key=tuple(np.round(z[:5],9))
                if key not in seen:poses.append(z);seen.add(key)
        z=seq[-1];key=tuple(np.round(z[:5],9))
        if key not in seen:poses.append(z);seen.add(key)
    eye=np.array([r['eye_local_ac_mm'] for r in existing['routes']]);ends=[]
    for z in poses:
        fk=v8.direct_kinematics(cfg,z[:5]);ends.append(eye@fk['body_rotation'].T+fk['ac_mm'])
    ends=np.array(ends)
    # All pair conflicts are cheap to evaluate in vectorized batches. Adding
    # them up front avoids repeatedly solving the same geometric alternatives.
    for g in range(8):
        ga=np.array([r['fixed_eye_world_mm'] for r in pools[g]])
        gb=ends[None,:,g,:]-5*(ends[None,:,g,:]-ga[:,None,:])/np.linalg.norm(ends[None,:,g,:]-ga[:,None,:],axis=2)[:,:,None]
        for h in range(g+1):
            hc=np.array([r['fixed_eye_world_mm'] for r in pools[h]])
            hd=ends[None,:,h,:]-5*(ends[None,:,h,:]-hc[:,None,:])/np.linalg.norm(ends[None,:,h,:]-hc[:,None,:],axis=2)[:,:,None]
            for offset in range(0,len(ga),4):
                dist=segment_distances(ga[offset:offset+4,None,None,:],gb[offset:offset+4,None,:,:],hc[None,:,None,:],hd[None,:,:,:]).min(2)
                for j,l in np.argwhere(dist<1.8):
                    a=int(starts[g]+offset+j);b=int(starts[h]+l)
                    if a>b:add([(a,1),(b,1)],-np.inf,1);cuts.add((a,b))
    print('PAIR_CONFLICTS',len(cuts),flush=True)
    for iteration in range(35):
        mat=coo_matrix((vals,(rows,cols)),shape=(len(lower),V)).tocsc()
        result=milp(objective,integrality=intg,bounds=bounds,constraints=LinearConstraint(mat,np.array(lower),np.array(upper)),options={'time_limit':60.,'mip_rel_gap':.01})
        if result.x is None:print('NO_SOLUTION',iteration,result.message,flush=True);return
        chosen=[(g,j) for g in range(8) for j in range(starts[g],starts[g+1]) if result.x[j]>.5];bad=[]
        for ig,(g,j) in enumerate(chosen):
            for h,l in chosen[:ig]:
                a=np.broadcast_to(pools[g][j-starts[g]]['fixed_eye_world_mm'],(len(poses),3));c=np.broadcast_to(pools[h][l-starts[h]]['fixed_eye_world_mm'],(len(poses),3));b=ends[:,g].copy();dd=ends[:,h].copy()
                b-=5*(b-a)/np.linalg.norm(b-a,axis=1)[:,None];dd-=5*(dd-c)/np.linalg.norm(dd-c,axis=1)[:,None]
                gap=float(segment_distances(a,b,c,dd).min())
                if gap<1.8:bad.append((j,l))
        print('MILP',iteration,'bad_pairs',len(bad),'solver',result.message,flush=True)
        if not bad:
            output={'pass':True,'selection_method':'mixed_integer selection; positive closure at 9 representative poses; physical route pools checked at 84 poses','pool_sizes':sizes,'terminal_clamp_length_mm':5.,'routes':[],'dense_drive_and_geometry_pass':False}
            attached=set();extra=9
            for g,j in chosen:
                rid=f'S{g+1:02d}' if g not in attached else f'S{extra:02d}'
                if g in attached:extra+=1
                attached.add(g)
                output['routes'].append({**existing['routes'][g],**pools[g][j-starts[g]],'id':rid,'attachment_id':f'S{g+1:02d}'})
            output['shoulder_wire_count']=len(chosen)
            (EXPORT/'full_arm_shoulder_routes_candidate.json').write_text(json.dumps(output,indent=2));print('FOUND',flush=True);return
        for a,b in bad:
            if (a,b) not in cuts:add([(a,1),(b,1)],-np.inf,1);cuts.add((a,b))
if __name__=='__main__':run()
