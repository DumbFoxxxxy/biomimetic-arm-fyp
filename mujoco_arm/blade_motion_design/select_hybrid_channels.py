"""Keep compact local channels where useful, retain clear external routes elsewhere."""
import copy,json,sys
import numpy as np
from scipy.optimize import linprog
from shoulder_model import SOURCE,EXPORT,v8
sys.stdout.reconfigure(encoding='utf-8')

def select():
    cfg=v8.load_config();fit=json.loads((SOURCE/'paper_motion_fit_repaired.json').read_text())
    base=json.loads((SOURCE/'anchor_search_no_contact_repaired.json').read_text())
    compact=json.loads((SOURCE/'anchor_search_no_contact_repaired_channels.json').read_text())
    # S07 was a genuinely local notch in the full sweep; other all-channel
    # corridors severed major regions and are not adopted blindly.
    result=copy.deepcopy(base); result['routes'][6]=copy.deepcopy(compact['routes'][6]);checks=[];eps=1e-6
    for pose in fit['poses']:
        q=np.array(pose['q_rad']);J=np.zeros((8,5))
        for k in range(5):
            p=q.copy();m=q.copy();p[k]+=eps;m[k]-=eps
            fp=v8.direct_kinematics(cfg,p);fm=v8.direct_kinematics(cfg,m)
            for j,r in enumerate(result['routes']):
                b=np.array(r['eye_local_ac_mm']);a=np.array(r['fixed_eye_world_mm'])
                J[j,k]=(np.linalg.norm(fp['ac_mm']+fp['body_rotation']@b-a)-np.linalg.norm(fm['ac_mm']+fm['body_rotation']@b-a))/(2*eps)
        W=J.T;W=W/np.linalg.norm(W,axis=0)
        lp=linprog(np.zeros(8),A_eq=np.vstack([W,np.ones((1,8))]),b_eq=np.r_[np.zeros(5),1.],bounds=[(1e-5,None)]*8,method='highs')
        checks.append({'id':pose['id'],'closure':bool(lp.success),'rank':int(np.linalg.matrix_rank(W))})
    result.update({'scope':'Hybrid: S07 local free-clearance channel; other routes external; no assumed contact reactions',
        'pass':all(c['closure'] and c['rank']==5 for c in checks),'pose_checks':checks,'channel_cables':['S07'],
        'frozen':False,'requires_scapular_channels':True,'channel_geometry_verified':False})
    (SOURCE/'anchor_search_hybrid_channels.json').write_text(json.dumps(result,indent=2))
    print({'pass':result['pass'],'passed':sum(c['closure'] for c in checks),'total':len(checks)})
if __name__=='__main__':select()
