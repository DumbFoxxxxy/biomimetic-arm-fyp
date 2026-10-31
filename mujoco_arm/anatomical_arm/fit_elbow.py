"""Bounded engineering calibration of elbow axis against articular clearance.
This is not a clinical axis measurement. Leaves all source meshes unchanged.
"""
from build import *
from scipy.optimize import differential_evolution
from scipy.spatial.transform import Rotation

def run():
    bones,r=audit_source();E=np.array(r['landmarks_mm']['E']);axis=np.array(r['landmarks_mm']['elbow_axis']);h=bones['humerus'];fore=[bones['ulna'],bones['radius']];counter=0
    def params(x):return E+np.array(x[:3]),unit(axis+[0,x[3],x[4]])
    angles=np.radians(np.arange(10,146,15).tolist()+[145])
    def score(x,report=False):
        nonlocal counter
        c,a=params(x);vol=[]
        for angle in angles:
            R=Rotation.from_rotvec(a*angle).as_matrix();T=np.c_[R,c-R@c]
            vol.append(sum(max(0.,float((h^s.transform(T)).volume())) for s in fore))
        value=max(vol)+.0002*float(np.dot(x,x))
        counter+=1
        if counter%100==0:print('axis fit',counter,'objective',value,flush=True)
        return vol if report else value
    initial=score(np.zeros(5),True);sol=differential_evolution(score,[(-8,8)]*3+[(-.12,.12)]*2,seed=17,popsize=5,maxiter=16,tol=.001,polish=False)
    centre,ax=params(sol.x);data={'method':'bounded collision-informed engineering calibration; not an anatomical ground-truth fit','initial_E_mm':E.tolist(),'initial_axis':axis.tolist(),'E_mm':centre.tolist(),'axis':ax.tolist(),'displacement_mm':float(np.linalg.norm(centre-E)),'axis_change_deg':float(np.degrees(np.arccos(np.clip(ax@axis,-1,1)))),'sample_angles_deg':np.degrees(angles).tolist(),'initial_overlap_volumes_mm3':initial,'optimized_overlap_volumes_mm3':score(sol.x,True),'candidate_parameters':sol.x.tolist(),'requires_full_1deg_reaudit':True}
    dump(OUT/'elbow_calibration.json',data);print(data,flush=True)

if __name__=='__main__':run()
