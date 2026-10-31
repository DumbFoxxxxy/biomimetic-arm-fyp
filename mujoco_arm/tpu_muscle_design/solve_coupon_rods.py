"""Normalized discrete-rod compression assay, not a calibrated TPU FEM model.

Stretching and bending energies produce equilibrium node positions. The end
plates prescribe assay displacement, not a hand-imposed interior deformation.
E is factored out; no force in newtons or safe material strain is asserted.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.optimize import minimize
OUT=Path(__file__).resolve().parents[1]/'exports/tpu_muscle_design'
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tmp/shoulder_motion_design_20261029'))
from search_paper_anchors import segment_distances

def topology(cells=5,subdivisions=4):
    nodes=[];lookup={};branches=[]
    def node(p,shared=True):
        key=tuple(p)
        if shared and key in lookup:return lookup[key]
        i=len(nodes);nodes.append(np.array(p,dtype=float))
        if shared:lookup[key]=i
        return i
    for c in range(cells):
        centre=np.array([(c+.5)*15,0.,0.]);ci=node(centre)
        for x in [c*15.,(c+1)*15.]:
            for y in [-5.,5.]:
                for z in [-15.,15.]:
                    corner=np.array([x,y,z]);end=node(corner);ids=[ci]
                    for t in np.linspace(0,1,subdivisions+1)[1:-1]:ids.append(node(centre*(1-t)+corner*t,False))
                    branches.append(ids+[end])
    return np.array(nodes),np.array(branches)

def run():
    rest,branches=topology();edges=np.array([(a,b) for line in branches for a,b in zip(line[:-1],line[1:])]);a,b=edges.T;rest_length=np.linalg.norm(rest[b]-rest[a],axis=1);edge_dir=(rest[b]-rest[a])/rest_length[:,None]
    pairs=[]
    for k in range(len(branches)):
        for j in range(3):pairs.append((k*4+j,k*4+j+1,1.,1.))
    ends={}
    for k,line in enumerate(branches):ends.setdefault(line[0],[]).append((k*4,1.));ends.setdefault(line[-1],[]).append((k*4+3,-1.))
    for incident in ends.values():
        for i,(e,s) in enumerate(incident):
            for f,t in incident[:i]:pairs.append((e,f,s,t))
    pairs=np.array(pairs);ei=pairs[:,0].astype(int);ej=pairs[:,1].astype(int);sign=pairs[:,2]*pairs[:,3];cos0=np.sum(edge_dir[ei]*edge_dir[ej],axis=1)*sign
    internal=np.arange(len(pairs))<len(branches)*3
    clamp_edges=np.array([k*4+3 for k,line in enumerate(branches) if np.isclose(rest[line[-1],0],0) or np.isclose(rest[line[-1],0],75)])
    area=np.pi;inertia=np.pi/4;kb=inertia/((rest_length[ei]+rest_length[ej])*.5)
    fixed=np.zeros(rest.shape,dtype=bool);plate=(np.isclose(rest[:,0],0)|np.isclose(rest[:,0],75));fixed[plate]=True
    guides=np.array([line[0] for line in branches]);fixed[guides,1:]=True;free=~fixed
    def energy(points):
        delta=points[b]-points[a];length=np.linalg.norm(delta,axis=1);t=delta/length[:,None];stretch=length-rest_length
        value=float(np.sum(area/(2*rest_length)*stretch**2));de=area/rest_length*stretch;edgegrad=de[:,None]*t
        cosine=np.sum(t[ei]*t[ej],axis=1)*sign;error=cosine-cos0
        value+=float(np.sum(np.where(internal,kb*(1-cosine),kb/2*error**2)))
        dc=np.where(internal,-kb,kb*error)*sign
        dti=dc[:,None]*t[ej];dtj=dc[:,None]*t[ei]
        for indices,g in [(ei,dti),(ej,dtj)]:
            transformed=(g-t[indices]*np.sum(t[indices]*g,axis=1)[:,None])/length[indices,None]
            np.add.at(edgegrad,indices,transformed)
        diff=t[clamp_edges]-edge_dir[clamp_edges];kc=inertia/rest_length[clamp_edges];value+=float(np.sum(kc[:,None]*diff**2)/2)
        cg=kc[:,None]*diff;cg=(cg-t[clamp_edges]*np.sum(t[clamp_edges]*cg,axis=1)[:,None])/length[clamp_edges,None];np.add.at(edgegrad,clamp_edges,cg)
        grad=np.zeros_like(points);np.add.at(grad,a,-edgegrad);np.add.at(grad,b,edgegrad)
        return value,grad
    # An independent directional derivative catches wrong assembled gradients.
    rng=np.random.default_rng(42);test=rest+rng.normal(0,.03,rest.shape);direction=rng.normal(size=rest.shape);direction/=np.linalg.norm(direction);h=1e-5
    numerical=(energy(test+h*direction)[0]-energy(test-h*direction)[0])/(2*h);analytic=float(np.sum(energy(test)[1]*direction));gradient_error=abs(numerical-analytic)
    if gradient_error>1e-6:raise RuntimeError(('gradient mismatch',numerical,analytic))
    state=rest.copy();records=[];states=[]
    for strain in np.linspace(0,.50,21):
        boundary=rest.copy();boundary[np.isclose(rest[:,0],75),0]=75*(1-strain)
        initial=state.copy();initial[plate]=boundary[plate]
        def objective(vector):
            pts=boundary.copy();pts[free]=vector;v,g=energy(pts);return v,g[free]
        result=minimize(objective,initial[free],jac=True,method='L-BFGS-B',options={'maxiter':2500,'ftol':1e-12,'gtol':1e-7,'maxls':40})
        state=boundary.copy();state[free]=result.x;length=np.linalg.norm(state[b]-state[a],axis=1);gap=1e9;worst=None
        # Trim original shared-junction fusion zones from branch ends. Report
        # this exclusion explicitly; it is not a full self-contact certificate.
        trimmed=[]
        for line in branches:
            points=state[line].copy();left=4.;right=4.;segments=[]
            for p,q in zip(points[:-1],points[1:]):
                dl=np.linalg.norm(q-p)
                if left>=dl:left-=dl;continue
                pp=p+left/dl*(q-p);left=0.;segments.append((pp,q))
            reverse=[]
            for p,q in segments[::-1]:
                dl=np.linalg.norm(q-p)
                if right>=dl:right-=dl;continue
                qq=q-right/dl*(q-p);right=0.;reverse.append((p,qq))
            trimmed.append(reverse)
        for i,ss in enumerate(trimmed):
            for j,tt in enumerate(trimmed[:i]):
                for p,q in ss:
                    for u,v in tt:
                        distance=float(segment_distances(p[None,:],q[None,:],u[None,:],v[None,:])[0])-2.
                        if distance<gap:gap=distance;worst=[j,i]
        row={'cell_span_shortening_fraction':float(strain),'cell_span_mm':75*(1-float(strain)),'normalized_energy_E_factored_out':float(result.fun),'optimizer_success':bool(result.success),'free_gradient_inf':float(np.max(abs(result.jac))),'stationary_equilibrium_verified':bool(np.max(abs(result.jac))<1e-6),
             'maximum_axial_segment_strain':float(np.max(abs(length/rest_length-1))), 'branch_centreline_envelope_mm':(state.max(0)-state.min(0)).tolist(),'minimum_nonjunction_branch_surface_gap_mm':gap,'closest_branches':worst}
        records.append(row);states.append(state.copy());print(json.dumps(row),flush=True)
    OUT.mkdir(exist_ok=True,parents=True);np.savez_compressed(OUT/'coupon_rod_equilibria.npz',rest=rest,branches=branches,states=np.array(states),shortening=np.linspace(0,.5,21))
    report={'scope':'Displacement-controlled normalized discrete rods; no calibrated TPU law, outer surface, friction or contact forces. Negative gaps invalidate the noncontact equilibrium beyond first contact. Shared-junction 4 mm regions are excluded from the branch gap screen. Ring stiffness and joint torsion are unresolved. Not a verified muscle actuator or safe contraction limit.',
        'gradient_directional_check_absolute_error':gradient_error,'nodes':len(rest),'segments':len(edges),'first_sampled_nonjunction_contact_fraction':next((r['cell_span_shortening_fraction'] for r in records if r['minimum_nonjunction_branch_surface_gap_mm']<0),None),'rows':records,'tpu_actuation_verified':False}
    (OUT/'coupon_rod_compression_report.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':run()
