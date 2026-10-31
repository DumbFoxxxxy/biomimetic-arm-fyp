"""Full TPU surface-graph rigidity fit under a hypothetical guide/cap pose.

This is geometric compatibility optimization, NOT TPU finite elements or a
calibrated elastic/contact equilibrium. Guide poses are prescribed hypotheses,
not physical external anchors. No manufacturing source mesh is overwritten.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from project_runtime import WorkspacePath as Path

class RigidityFit:
    def __init__(self,V,faces,fixed):
        self.V=V;self.N=len(V);self.edges=np.unique(np.sort(np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);self.i,self.j=self.edges.T;self.dp=V[self.i]-V[self.j];self.w=1/np.maximum(np.linalg.norm(self.dp,axis=1),.2);self.fixed=np.flatnonzero(fixed);self.free=np.flatnonzero(~fixed);diagonal=np.bincount(np.r_[self.i,self.j],weights=np.r_[self.w,self.w],minlength=self.N);L=coo_matrix((np.r_[-self.w,-self.w,diagonal],(np.r_[self.i,self.j,np.arange(self.N)],np.r_[self.j,self.i,np.arange(self.N)])),shape=(self.N,self.N)).tocsc();self.Lfc=L[self.free][:,self.fixed];print('factor graph',self.N,'vertices',len(self.edges),'edges',len(self.fixed),'fixed',flush=True);self.lu=splu(L[self.free][:,self.free]);print('factor complete',flush=True)
    def rotations(self,X):
        difference=X[self.i]-X[self.j];cov=self.w[:,None,None]*difference[:,:,None]*self.dp[:,None,:];S=np.zeros((self.N,3,3));np.add.at(S,self.i,cov);np.add.at(S,self.j,cov);U,_,Vt=np.linalg.svd(S);bad=np.linalg.det(U@Vt)<0;U[bad,:,2]*=-1;return U@Vt
    def step(self,X,target):
        R=self.rotations(X);term=.5*self.w[:,None]*np.einsum('eij,ej->ei',R[self.i]+R[self.j],self.dp);b=np.zeros((self.N,3));np.add.at(b,self.i,term);np.add.at(b,self.j,-term);Y=np.empty_like(X);Y[self.fixed]=target[self.fixed];Y[self.free]=self.lu.solve(b[self.free]-self.Lfc@Y[self.fixed]);return Y,R
    def energy(self,X,R):
        delta=X[self.i]-X[self.j];a=delta-np.einsum('eij,ej->ei',R[self.i],self.dp);b=delta-np.einsum('eij,ej->ei',R[self.j],self.dp);return float((.5*self.w*((a*a).sum(1)+(b*b).sum(1))).sum())

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';candidate=parent/'captive_keepers';study_path=candidate/'wire_pose_regularized_fine/wire_guides_single_132.json';study=json.loads(study_path.read_text());case=study['poses'][0];assert case['muscle']=='FJ1486' and case['elbow_deg']==132
    for p,h in study['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    a=json.loads((candidate/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};fid=case['muscle'];folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());r=rows[fid+'_TPU_muscle'];solid,V,faces=load(ROOT/r['file']);V=V+np.array(r['translation_display_mm']);offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));A,B=np.array(d['pad_centres_display_mm'])+offset;F=np.array(d['basis_columns_display']);z=F[:,2];ha=(V-A)@F;hb=(V-B)@F
    # Dense pads and material immediately around the four existing M3 clamp
    # seats follow actual caps; central soft transitions are not fixed in space.
    near_bolt_a=np.minimum((ha[:,0]-7)**2+ha[:,1]**2,(ha[:,0]+7)**2+ha[:,1]**2)<4**2;near_bolt_b=np.minimum((hb[:,0]-7)**2+hb[:,1]**2,(hb[:,0]+7)**2+hb[:,1]**2)<4**2;prox=(ha[:,2]<=2)|(near_bolt_a&(ha[:,2]<=5.4));dist=(hb[:,2]>=-2)|(near_bolt_b&(hb[:,2]>=-5.4));fixed=prox|dist;guide_masks=[]
    for g in d['actual_short_tube_guides']:
        c=np.array(g['centre_display_mm'])+offset;axis=np.array(g['axis_display']);delta=V-c;axial=delta@axis;radial2=np.maximum(0,(delta*delta).sum(1)-axial*axial);mask=(abs(axial)<=g['length_mm']/2+1e-7)&(radial2<=(g['outer_radius_mm']+.015)**2);assert not np.any(mask&(prox|dist));fixed|=mask;guide_masks.append(mask)
    dest=candidate/'lattice_graph_probe';dest.mkdir(exist_ok=True);inputs=[study_path,candidate/'assembly.json',folder/'design.json',ROOT/r['file'],parent/'conditional_motion.json']+[ROOT/item['file'] for item in rows.values() if item['kind']!='TPU'];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};solver=RigidityFit(V,faces,fixed)
    # Independent mathematical sanity gates on the actual whole manufacturing
    # mesh: exact rest and a known Euclidean rigid transform must be preserved.
    rest,_=solver.step(V,V);rest_error=float(np.linalg.norm(rest-V,axis=1).max());assert rest_error<1e-6,rest_error
    from scipy.spatial.transform import Rotation
    Rtest=Rotation.from_rotvec(np.array([.1,-.2,.3])).as_matrix();rigid=V@Rtest.T+[3,-4,6];recovered,_=solver.step(rigid,rigid);rigid_error=float(np.linalg.norm(recovered-rigid,axis=1).max());assert rigid_error<1e-6,rigid_error
    pose=next(p for p in json.loads((parent/'conditional_motion.json').read_text())['poses'] if p['elbow_deg']==132);Tc=np.array(pose['cap_poses'][fid]);Tf=np.array(pose['forearm_pose']);target=V.copy();target[dist]=V[dist]@Tc[:,:3].T+Tc[:,3]
    for mask,g in zip(guide_masks,case['guide_checks']):T=np.array(g['pose_affine_mm']);target[mask]=V[mask]@T[:,:3].T+T[:,3]
    u=np.clip(((V-A)@z)/float((B-A)@z),0,1);w=6*u**5-15*u**4+10*u**3;X=V+w[:,None]*(V@Tc[:,:3].T+Tc[:,3]-V);X[fixed]=target[fixed];history=[];converged=False
    for iteration in range(80):
        Y,R=solver.step(X,target);change=float(np.linalg.norm(Y-X,axis=1).max());energy=solver.energy(Y,R);history.append({'iteration':iteration+1,'maximum_vertex_change_mm':change,'geometric_energy_nonphysical_units':energy});X=Y
        if iteration%10==0:print('rigidity iteration',iteration+1,'maxchange',change,'energy',energy,flush=True)
        if change<.005:converged=True;break
    final_R=solver.rotations(X);edge_new=np.linalg.norm(X[solver.i]-X[solver.j],axis=1);edge_old=np.linalg.norm(solver.dp,axis=1);valid=edge_old>.05;ratio=edge_new[valid]/edge_old[valid];mesh=mf.Manifold(mf.Mesh64(X.astype(np.float64),faces.astype(np.uint64)));contacts=[];wire_hit=None
    if mesh.status()==mf.Error.NoError:
        for key,item in rows.items():
            if item['kind']=='TPU':continue
            other=load(ROOT/item['file'])[0].translate(item['translation_display_mm']);owner=item['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;other=other.transform(T) if T is not None else other;v=intersection(mesh,other)
            if v>.001:contacts.append({'part':key,'volume_mm3':v,'classification':'geometric_hypothesis_contact_requires_solution'})
        from study_swivel_wire_and_guides import capsule_chain
        wire_hit=intersection(mesh,capsule_chain(np.array(case['centreline_display_mm']),.3))
    target_error=float(np.linalg.norm(X[fixed]-target[fixed],axis=1).max());report={'scope':'whole manufactured TPU surface-graph rigidity fit with hypothetical guide poses andactual cap clamp regions; not3DsolidFEM orcalibratedTPU equilibrium/contact/friction proof','muscle':fid,'elbow_deg':132,'vertices':len(V),'triangles':len(faces),'surface_graph_edges':len(solver.edges),'fixed_material_vertices':int(fixed.sum()),'rest_patch_error_mm':rest_error,'rigid_transform_patch_error_mm':rigid_error,'iterations':history,'iteration_tolerance_converged':converged,'material_constraints_error_mm':target_error,'surface_edge_stretch_min_max':[float(ratio.min()),float(ratio.max())],'surface_edge_stretch_p01_p99':np.percentile(ratio,[1,99]).tolist(),'mapped_mesh_status':str(mesh.status()),'material_volume_change_fraction':float(mesh.volume()/solid.volume()-1) if mesh.status()==mf.Error.NoError else None,'rigid_contact_reviews':contacts,'steel_TPU_intersection_mm3':wire_hit,'full_self_contact_verified':False,'actual_physical_TPU_or_driverange_verified':False,'source_sha256':hashes,'frozen':False};(dest/'graph_review.json').write_text(json.dumps(report,indent=2),encoding='utf8');(dest/'TPU_GEOMETRIC_HYPOTHESIS_NOT_PRINT.obj').write_text('# GEOMETRIC COMPATIBILITY HYPOTHESIS ONLY. Not a manufacturing part or calibrated TPU simulation.\n'+'\n'.join(['v '+' '.join(format(float(c),'.17g') for c in p) for p in X]+['f '+' '.join(str(int(i)+1) for i in face) for face in faces])+'\n',encoding='utf8');assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());print('graph probe complete contacts',len(contacts),'wire',wire_hit,'volumechange',report['material_volume_change_fraction'],flush=True)

if __name__=='__main__':run()
