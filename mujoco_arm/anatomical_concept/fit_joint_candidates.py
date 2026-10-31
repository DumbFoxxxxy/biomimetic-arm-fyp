"""Fit engineering joint candidates from current raw anatomy and support geometry.
Not imported from V2/A1. No anatomical surface is altered or scaled.
"""
from build_reference import *
from extract_attachment_regions import closest
from scipy.optimize import least_squares,differential_evolution
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(ROOT/'tmp/shoulder_motion_design_20261029/vendor'))
import manifold3d as mf

def unit(x):return np.asarray(x)/np.linalg.norm(x)

def sphere(points):
    x=np.linalg.lstsq(np.c_[2*points,np.ones(len(points))],(points*points).sum(1),rcond=None)[0];r=np.sqrt(max(0,x[3]+x[:3]@x[:3]));fit=least_squares(lambda a:np.linalg.norm(points-a[:3],axis=1)-a[3],np.r_[x[:3],r],loss='soft_l1');res=np.linalg.norm(points-fit.x[:3],axis=1)-fit.x[3];return fit.x[:3],float(fit.x[3]),float(np.sqrt(np.mean(res*res)))

def run():
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']};meshes={}
    for fid in ['FJ3362','FJ3384','FJ3368','FJ3391','FJ3349','FJ3290']:
        v,f=read_obj(str(ROOT/byid[fid]['source']));v,f=weld_stl(v[f]);meshes[fid]=(v,f)
    clav,_=meshes['FJ3362'];scap,sf=meshes['FJ3384'];sternum,tf=meshes['FJ3290'];h,hf=meshes['FJ3368'];ulna,uf=meshes['FJ3391'];radius,rf=meshes['FJ3349']
    medial=clav[clav[:,0]>=np.percentile(clav[:,0],96)].mean(0);lateral=clav[clav[:,0]<=np.percentile(clav[:,0],4)].mean(0);sq,sd,_=closest(medial,sternum[tf]);aq,ad,_=closest(lateral,scap[sf]);SC=(medial+sq)/2;AC=(lateral+aq)/2
    mask=(h[:,2]>h[:,2].max()-36)&(h[:,0]>h[:,0].max()-38);GH,r,rms=sphere(h[mask]);variants=[]
    for zband,xband in [(32,34),(36,38),(40,42)]:
        p=h[(h[:,2]>h[:,2].max()-zband)&(h[:,0]>h[:,0].max()-xband)];c,rr,err=sphere(p);variants.append({'z_band_mm':zband,'medial_x_band_mm':xband,'centre_source_mm':c.tolist(),'radius_mm':rr,'radial_rms_mm':err,'centre_shift_from_selected_mm':float(np.linalg.norm(c-GH))})
    distal=h[h[:,2]<h[:,2].min()+25];E0=distal.mean(0);_,_,vt=np.linalg.svd(distal-E0,full_matrices=False);axis0=vt[0];axis0*=1 if axis0[0]<0 else -1
    H=mf.Manifold(mf.Mesh64(h.astype(np.float64),hf.astype(np.uint64)));fore=[mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64))) for v,f in [(ulna,uf),(radius,rf)]];angles=np.radians(list(range(10,151,10)));calls=0
    def params(x):return E0+np.array(x[:3]),unit(axis0+[0,x[3],x[4]])
    def score(x,detail=False):
        nonlocal calls
        E,axis=params(x);vol=[]
        for a in angles:
            R=Rotation.from_rotvec(axis*a).as_matrix();T=np.c_[R,E-R@E];vol.append(sum(max(0,float((H^s.transform(T)).volume())) for s in fore))
        calls+=1
        if calls%150==0:print('joint fit',calls,flush=True)
        return vol if detail else max(vol)+.0002*float(np.dot(x,x))
    initial=score(np.zeros(5),True);fit=differential_evolution(score,[(-8,8)]*3+[(-.12,.12)]*2,popsize=5,maxiter=18,seed=29,polish=False);E,axis=params(fit.x);sweep=[];first=None
    for a in range(151):
        R=Rotation.from_rotvec(axis*np.radians(a)).as_matrix();T=np.c_[R,E-R@E];vol=sum(max(0,float((H^s.transform(T)).volume())) for s in fore);sweep.append({'angle_deg':a,'overlap_mm3':vol})
        if vol>.001 and first is None:first=a
    data={'scope':'geometry-derived engineering candidates, not clinical landmarks or installed joints','coordinate_system':'raw BodyParts3D mm, no old frames or offsets','source_hashes':{fid:byid[fid]['sha256'] for fid in meshes},'SC':{'centre_source_mm':SC.tolist(),'clavicle_medial_end_mean_mm':medial.tolist(),'manubrium_nearest_surface_mm':sq.tolist(),'sample_end_gap_mm':sd,'method':'midpoint between clavicle medial4% vertex-band centroid and manubrium surface'},'AC':{'centre_source_mm':AC.tolist(),'clavicle_lateral_end_mean_mm':lateral.tolist(),'scapula_nearest_surface_mm':aq.tolist(),'sample_end_gap_mm':ad,'method':'midpoint between lateral4% clavicle-band centroid and scapula surface'},'GH':{'centre_source_mm':GH.tolist(),'radius_mm':r,'radial_rms_mm':rms,'selection':'proximal36mm and medial38mm humerus bands','selection_sensitivity':variants},'elbow':{'initial_centre_source_mm':E0.tolist(),'initial_axis_source':axis0.tolist(),'centre_source_mm':E.tolist(),'axis_source':axis.tolist(),'centre_adjustment_mm':float(np.linalg.norm(E-E0)),'axis_adjustment_deg':float(np.degrees(np.arccos(np.clip(axis@axis0,-1,1)))),'initial_overlap_samples_mm3':initial,'calibration_overlap_samples_mm3':score(fit.x,True),'checked_bone_sweep':sweep,'checked_until_deg':150,'first_colliding_integer_sample_deg':first,'volume_threshold_mm3':.001,'forearm_assumption':'ulna/radius held in source mutual pose for this elbow-only fit; this does not decide final pronation/supination design'},'proposed_constraint_candidate':'SC spherical3 + AC spherical3 + GH spherical3 + elbow hinge1. Need physical retention/support and muscle controllability checks; not frozen.','not_verified':['clinical landmark accuracy','SC/AC/GH motion collisions','hardware interference and retention','motor/tendon controllability','TPU deformation and force'],'frozen':False};(OUT/'joint_candidates.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf8');print('SC/AC surface gaps',sd,ad,'GH radius/rms',r,rms,'elbow adjust',data['elbow']['centre_adjustment_mm'],'first collision',first,flush=True)

if __name__=='__main__':run()
