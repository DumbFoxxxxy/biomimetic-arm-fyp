"""Search GH cuff attachments over the full requested single-axis motions.

No payload criterion. Straight candidate routes must clear the joint-center
keepout; full bone/cable geometry is a separate verification stage.
"""
import json,sys,itertools,argparse
import numpy as np
from scipy.optimize import linprog
from shoulder_model import EXPORT,v8,SOURCE
sys.path.insert(0,str(SOURCE))
from search_paper_anchors import segment_distances,closest_surface_point
sys.stdout.reconfigure(encoding='utf-8')
A=np.array([.087155743,-.996194698,0.]);A/=np.linalg.norm(A)
B=np.array([.996194698,.087155743,0.]);B/=np.linalg.norm(B)
C=np.array([0.,0.,1.])

def rotation(q):
    return v8.rotation_matrix(B,q[0])@v8.rotation_matrix(C,q[1])@v8.rotation_matrix(A,q[2])

def poses():
    rows=[{'id':'ZERO','q_deg':[0.,0.,0.]}]
    for name,k,end in [('FLEX',0,120.),('EXT',0,-40.),('ABD',1,120.),('INT',2,70.),('EXT_ROT',2,-90.)]:
        for f in np.linspace(.05,1.,20):
            q=np.zeros(3);q[k]=end*f;rows.append({'id':f'{name}_{f:.2f}','q_deg':q.tolist()})
    return rows

def search(geometry=False):
    target=poses();Rs=np.array([rotation(np.radians(p['q_deg'])) for p in target]);n=len(Rs);groups=[]
    root_tri=None
    if geometry:
        mesh=np.load(EXPORT/'scapula_hybrid_design_geometry.npz');gh_from_ac=np.array([159.06,187.54,7.30])-np.array(v8.load_config()['geometry_mm']['AC'])
        root_tri=mesh['vertices_ac_local_mm'][mesh['faces']]-gh_from_ac
    shaft0=np.einsum('nij,j->ni',Rs,15*A);shaft1=np.einsum('nij,j->ni',Rs,305*A)
    for g in range(8):
        phi=2*np.pi*g/8;b=(35. if geometry else 25.)*(np.cos(phi)*B+np.sin(phi)*C)+45.*A
        moving=np.einsum('nij,j->ni',Rs,b);pool=[]
        for psi,height,radius in itertools.product(np.linspace(0,2*np.pi,24,endpoint=False),[-60.,-30.,0.,30.],[40.,60.,80.]):
            a=radius*(np.cos(psi)*B+np.sin(psi)*C)+height*A
            delta=moving-a;tt=np.clip(-np.sum(a*delta,axis=1)/np.sum(delta*delta,axis=1),0,1)
            clear=float(np.min(np.linalg.norm(a+tt[:,None]*delta,axis=1)))
            if clear<12.:continue
            support=None
            if geometry:
                unit=(moving-a)/np.linalg.norm(moving-a,axis=1)[:,None]
                free_end=moving-3.*unit
                if np.min(segment_distances(np.broadcast_to(a,(n,3)),free_end,shaft0,shaft1))<21.:continue
                support=closest_surface_point(a,root_tri)
                if np.min(segment_distances(np.broadcast_to(support,(n,3)),np.broadcast_to(a,(n,3)),shaft0,shaft1))<22.:continue
                delta=a-support;tt=np.clip(-np.dot(support,delta)/max(np.dot(delta,delta),1e-10),0,1)
                if np.linalg.norm(support+tt*delta)<12.:continue
            length=np.linalg.norm(a-moving,axis=1);direction=(a-moving)/length[:,None]
            W=np.cross(moving,direction);norm=np.linalg.norm(W,axis=1)
            if np.any(norm<1e-4):continue
            pool.append({'origin_local_gh_mm':a.tolist(),'insertion_local_gh_mm':b.tolist(),'W':W/norm[:,None],
                'joint_keepout_distance_mm':clear,'length':length,
                'support_root_local_gh_mm':None if support is None else support.tolist()})
        groups.append(pool);print('GH_GROUP',g+1,'candidates',len(pool),flush=True)
    if any(not p for p in groups):
        report={'pass':False,'failure':'empty_straight_route_pool','pool_counts':[len(p) for p in groups]}
    else:
        all_pool=[r for group in groups for r in group]
        coverage=[]
        for index in range(n):
            W=np.column_stack([r['W'][index] for r in all_pool])
            lp=linprog(np.zeros(len(all_pool)),A_eq=np.vstack([W,np.ones((1,len(all_pool)))]),b_eq=[0,0,0,1],bounds=[(1e-7,None)]*len(all_pool),method='highs')
            coverage.append(bool(lp.success))
        print('ALL_CANDIDATE_COVERAGE',sum(coverage),'/',n,flush=True)
        if geometry:
            np.savez_compressed(EXPORT/'gh_geometry_pool.npz',W=np.array([r['W'] for r in all_pool]))
            (EXPORT/'gh_geometry_pool.json').write_text(json.dumps({'poses':target,'groups':[len(g) for g in groups],
                'routes':[{k:v for k,v in r.items() if k not in ['W','length']} for r in all_pool],
                'all_pool_coverage':[{'id':p['id'],'pass':okay} for p,okay in zip(target,coverage)]},indent=2))
        indices=[0,20,40,60,80,100];rng=np.random.default_rng(301);best=None;bestscore=-1;passing=None
        def check(selection,index):
            W=np.column_stack([groups[g][j]['W'][index] for g,j in enumerate(selection)])
            if np.linalg.matrix_rank(W)<3:return False
            lp=linprog(np.zeros(8),A_eq=np.vstack([W,np.ones((1,8))]),b_eq=[0,0,0,1],bounds=[(1e-5,None)]*8,method='highs')
            return bool(lp.success)
        for trial in range(6000):
            if best is not None and rng.random()<.8:
                selection=list(best)
                for g in rng.choice(8,size=int(rng.integers(1,4)),replace=False):selection[g]=int(rng.integers(len(groups[g])))
            else:selection=[int(rng.integers(len(p))) for p in groups]
            score=sum(check(selection,i) for i in indices)
            if score>bestscore or (score==bestscore and rng.random()<.03):best=selection;bestscore=score
            if score==len(indices) and all(check(selection,i) for i in range(n)):passing=selection;break
        selected=passing or best;routes=[]
        for g,j in enumerate(selected):
            row=groups[g][j];routes.append({'id':f'G{g+1:02d}',**{k:v for k,v in row.items() if k not in ['W','length']},
                'minimum_length_mm':float(row['length'].min()),'maximum_length_mm':float(row['length'].max()),
                'excursion_mm':float(np.ptp(row['length']))})
        report={'pass':passing is not None,'scope':'GH straight routes, positive torque spanning and 12mm center keepout only',
                'trial_count':trial+1,'routes':routes,'poses':target,'full_geometry_verified':False,
                'pass_count':sum(check(selected,i) for i in range(n)),'pose_count':n,'all_pool_pass_count':sum(coverage)}
    report['conservative_shaft_and_support_filter']=geometry
    EXPORT.mkdir(parents=True,exist_ok=True);(EXPORT/('gh_anchor_search_geometry.json' if geometry else 'gh_anchor_search.json')).write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['routes','poses']},indent=2),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--geometry',action='store_true');args=parser.parse_args();search(args.geometry)
