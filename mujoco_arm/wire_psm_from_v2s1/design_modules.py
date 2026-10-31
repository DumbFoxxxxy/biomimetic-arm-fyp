"""PSM belly on one carrier, mobile belly end C, passive tendon C-G-B.

The present search screens the maximum A-G-B steel corridor. Actual C
positions and TPU envelope are evaluated separately, never assumed identical
to wire take-up. Polished spherical R8 guides are geometric candidates.
"""
from common import *
from fit_mounts import outer_triangles
from audit_interfaces import wire_pair_gap
from scipy.optimize import linprog
import manifold3d as md
import csv

def run():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--only',default='');parser.add_argument('--compact',action='store_true');parser.add_argument('--no-pair-check',action='store_true');parser.add_argument('--preserve-closure',action='store_true');parser.add_argument('--space-check',action='store_true');args=parser.parse_args();only=set(args.only.split(',')) if args.only else None
    cached={};cached_root=None
    if only:
        cached={z['id']:z for z in json.loads((OUT/'psm_module_candidate.json').read_text())['modules']};cached_root=ET.parse(OUT/'psm_module_candidate.xml').getroot()
    rng=np.random.default_rng(20261029);oldtree=template();oldm,oldd=compile_tree(oldtree);oldroot=oldtree.getroot()
    tree=ET.parse(BASE/'skeleton.xml');root=tree.getroot();parts_model,parts_data=compile_tree(tree);parts=solids(parts_model,parts_data)
    byname={mujoco.mj_id2name(parts_model,mujoco.mjtObj.mjOBJ_BODY,p['body']):(n,p) for n,p in parts.items()}
    pool={}
    for body,(name,part) in byname.items():
        tri=outer_triangles(name,part);centers=tri.mean(1);normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(normals,axis=1);normals/=np.maximum(area,1e-20)[:,None]
        if args.compact:
            ids=rng.choice(len(tri),min(700,len(tri)),replace=False,p=area/area.sum());pool[body]=(centers[ids],normals[ids])
        else:
            ids=rng.choice(len(tri),1400,replace=True,p=area/area.sum());r=np.sqrt(rng.random(len(ids)));t=rng.random(len(ids));sample=(1-r)[:,None]*tri[ids,0]+(r*(1-t))[:,None]*tri[ids,1]+(r*t)[:,None]*tri[ids,2];pool[body]=(sample,normals[ids])
    tendon=ET.SubElement(root,'tendon');modules=[]
    for oldroute in oldroot.find('tendon'):
        rid=oldroute.get('name').replace('cable_','');oldsites=[n.get('site') for n in oldroute if n.tag=='site'];sid0=mujoco.mj_name2id(oldm,mujoco.mjtObj.mjOBJ_SITE,oldsites[0]);sid1=mujoco.mj_name2id(oldm,mujoco.mjtObj.mjOBJ_SITE,oldsites[-1]);origin=oldd.site_xpos[sid0]*1000;insertion=oldd.site_xpos[sid1]*1000
        if rid.startswith('E'):carrier='humerus';target='ulna';seed=oldd.xpos[oldmujocoid(oldm,'ulna')]*1000+np.array([3.92,15 if rid=='E01' else -15,-44.8]);gseed=oldd.xpos[oldmujocoid(oldm,'humerus')]*1000+np.array([24,20 if rid=='E01' else -20,-250])
        elif rid in ['G01','G02','G03','G04']:
            carrier='humerus';target='clavicle_body' if rid=='G01' else 'scapula';gh=oldd.xpos[oldmujocoid(oldm,'humerus')]*1000
            seed=(oldd.xpos[oldmujocoid(oldm,'scapula')]*1000+np.array([-15,10,0])) if rid=='G01' else gh+({'G02':np.array([-15,-20,-15]),'G03':np.array([0,-15,15]),'G04':np.array([-70,-55,-70])}[rid]);gseed=gh+np.array([4,0,-45])
        elif rid.startswith('G'):
            carrier='scapula';target='humerus';gh=oldd.xpos[oldmujocoid(oldm,'humerus')]*1000;seed=gh+np.array([4.35,15 if rid=='G05' else -15,-49.8]);gseed=gh+np.array([-12,-20,0])
        elif rid.startswith('S'):
            carrier='thorax';target='scapula';seed=insertion;gseed=origin
            if rid=='S01':gseed=np.array([75.,-95.,265.]);seed=np.array([70.,-100.,150.])
            elif rid=='S03':gseed=np.array([135.,-45.,230.]);seed=np.array([65.,-90.,190.])
        else:
            carrier='thorax';target='humerus';seed=insertion;guide=[n for n in oldroute if n.tag=='geom'];gseed=oldd.geom_xpos[mujoco.mj_name2id(oldm,mujoco.mjtObj.mjOBJ_GEOM,guide[0].get('geom'))]*1000
        for role,body in [('A',carrier),('C',carrier),('B',target),('Groot',carrier),('Aroot',carrier),('Broot',target)]:ET.SubElement(root.find(f'.//body[@name="{body}"]'),'site',{'name':rid+'_'+role,'pos':'0 0 0','size':'.002','rgba':'0.9 0.4 0.1 1'})
        ET.SubElement(root.find(f'.//body[@name="{carrier}"]'),'geom',{'name':rid+'_guide','type':'sphere','size':'.008','pos':'0 0 0','rgba':'.7 .7 .75 1','contype':'0','conaffinity':'0','density':'0'})
        path=ET.SubElement(tendon,'spatial',{'name':rid+'_passive','width':'.0005','rgba':'.85 .3 .1 1'});ET.SubElement(path,'site',{'site':rid+'_C'});ET.SubElement(path,'geom',{'geom':rid+'_guide'});ET.SubElement(path,'site',{'site':rid+'_B'})
        modules.append({'id':rid,'carrier':carrier,'insertion_body':target,'gseed_world':gseed,'bseed_world':seed})
    if only:
        for item in modules:
            rid=item['id']
            for role in ['A','C','B','Aroot','Groot','Broot']:
                root.find(f'.//site[@name="{rid}_{role}"]').set('pos',cached_root.find(f'.//site[@name="{rid}_{role}"]').get('pos'))
            root.find(f'.//geom[@name="{rid}_guide"]').set('pos',cached_root.find(f'.//geom[@name="{rid}_guide"]').get('pos'))
    m,d=compile_tree(tree);coarse=[np.zeros(7)]
    for seq in pose_library().values():
        chosen=np.r_[np.linspace(0,len(seq)-1,min(4,len(seq))).astype(int),[i for i in [1,2,3,5] if i<len(seq)]];coarse.extend(np.array(seq)[np.unique(chosen)])
    old_audit=OUT/'interfaces_geometry_audit.json'
    if old_audit.exists():
        for action in json.loads(old_audit.read_text()).get('actions',[]):
            for field in ['first_hard_geometry_failure','first_wire_pair_failure']:
                if action.get(field):coarse.append(np.array(action[field]['q_deg']))
    coarse=np.unique(np.array(coarse),axis=0);coarse=np.vstack([np.zeros(7),coarse[np.any(np.abs(coarse)>1e-12,axis=1)]])
    results=[]
    for index,module in enumerate(modules):
        if only and module['id'] not in only:
            rid=module['id']
            for role in ['A','C','B','Aroot','Groot','Broot']:
                root.find(f'.//site[@name="{rid}_{role}"]').set('pos',cached_root.find(f'.//site[@name="{rid}_{role}"]').get('pos'))
            root.find(f'.//geom[@name="{rid}_guide"]').set('pos',cached_root.find(f'.//geom[@name="{rid}_guide"]').get('pos'));results.append(cached[rid]);m,d=compile_tree(tree);continue
        rid=module['id'];carrier=module['carrier'];target=module['insertion_body'];cb=oldmujocoid(m,carrier);tb=oldmujocoid(m,target);pc,nc=pool[carrier];pt,nt=pool[target]
        seedg=module['gseed_world']-d.xpos[cb]*1000;seedb=module['bseed_world']-d.xpos[tb]*1000
        gi=np.argsort(np.linalg.norm(pc-seedg,axis=1))[:min(120,len(pc))];bi=np.argsort(np.linalg.norm(pt-seedb,axis=1))[:min(140,len(pt))]
        if rid in ['G05','G06']:
            wanted=np.array([0,1 if rid=='G05' else -1,0]);keep=bi[(nt[bi]@wanted)>.5]
            if len(keep):bi=keep
        if carrier=='humerus':
            wanted=np.array([0,1,0]) if rid in ['E01','G01'] else np.array([0,-1,0]) if rid in ['E02','G02'] else np.array([1,0,0]) if rid=='G03' else np.array([-1,0,0])
            keep=gi[(nc[gi]@wanted)>.35]
            if len(keep):gi=keep
        ids={role:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,rid+'_'+role) for role in ['A','C','B','Aroot','Groot','Broot']};gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,rid+'_guide');best=None
        target_axis=6 if rid.startswith('E') else 3 if rid in ['G01','G02','P01','P02'] else 4 if rid in ['G03','G04'] else 5 if rid in ['G05','G06'] else {'S01':0,'S02':0,'S03':1,'S04':1,'S05':2,'S06':2}[rid]
        target_sign=-1 if rid in ['E02','G02','G04','G06','S02','S04','S06'] else 1
        reference_poses=[]
        for reference_q in coarse:
            d.qpos[:]=np.radians(reference_q);mujoco.mj_forward(m,d);references=[]
            for j,other in enumerate(modules):
                if args.no_pair_check or j==index or (not only and j>=index):continue
                oid=other['id'];ai=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,oid+'_A');bi2=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,oid+'_B');gg=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,oid+'_guide');rr=swivel_route(d.site_xpos[ai]*1000,d.geom_xpos[gg]*1000,d.site_xpos[bi2]*1000)
                if rr is not None and rr['wrap_deg']<=180:references.append(rr['path'])
            reference_poses.append(references)
        preferred_AG=200 if carrier=='humerus' else 120 if carrier=='scapula' else 190
        old_matrix=None
        if args.preserve_closure:
            old_matrix=np.zeros((16,7))
            for joint in range(7):
                lengths=[]
                for sign in [1,-1]:
                    d.qpos[:]=0;d.qpos[joint]=sign*1e-6;mujoco.mj_forward(m,d);vals=[]
                    for other in modules:
                        oid=other['id'];ai=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,oid+'_A');bi2=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_SITE,oid+'_B');gg=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,oid+'_guide');rr=swivel_route(d.site_xpos[ai]*1000,d.geom_xpos[gg]*1000,d.site_xpos[bi2]*1000);vals.append(rr['D_mm'] if rr is not None else np.nan)
                    lengths.append(np.array(vals))
                old_matrix[:,joint]=-(lengths[0]-lengths[1])/2e-6
            old_matrix[np.abs(old_matrix)<1e-4]=0
        for trial in range(700 if args.compact else 500):
            gidx=gi[trial%len(gi)];bidx=bi[(trial*17)%len(bi)];gfoot=pc[gidx];bfoot=pt[bidx];gn=nc[gidx];bn=nt[bidx]
            length=np.linalg.norm(pc-gfoot,axis=1);align=nc@gn;choices=np.where((length>85)&(length<225)&(align>.25))[0]
            if not len(choices):continue
            aidx=choices[(trial*23)%len(choices)];afoot=pc[aidx];an=nc[aidx];offset=[20.,28.,36.,44.][trial%4];A=afoot+an*offset;G=gfoot+gn*offset;B=bfoot+bn*[8.,16.,24.][trial%3]
            # User steering: length is a preference, not a hard requirement.
            # Try longer compartments briefly, then allow actual short space.
            if np.linalg.norm(G-A)<(75 if args.compact else preferred_AG if trial<160 else 85):continue
            m.site_pos[ids['A']]=A*.001;m.site_pos[ids['C']]=A*.001;m.site_pos[ids['B']]=B*.001;m.geom_pos[gid]=G*.001
            def route_at(q):
                d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);return swivel_route(d.site_xpos[ids['A']]*1000,d.geom_xpos[gid]*1000,d.site_xpos[ids['B']]*1000)
            plus=np.zeros(7);plus[target_axis]=1e-4;minus=-plus;rp=route_at(plus);rm=route_at(minus)
            moment=float(-(rp['D_mm']-rm['D_mm'])/(2*np.radians(1e-4))) if rp is not None and rm is not None else -1e6;functional=bool(target_sign*moment>=8.)
            if rid in ['G05','G06']:functional=functional and target_sign*moment<=18
            clear=0;zero=False;worst=1e6
            for qi,q in enumerate(coarse):
                routed=route_at(q);blocked=routed is None or routed['wrap_deg']>180
                if blocked:break
                path=routed['path']
                for bone,part in parts.items():
                    local=(path-d.xpos[part['body']]*1000)@d.xmat[part['body']].reshape(3,3);gap=local_gap(local,part);worst=min(worst,gap)
                    if gap<(.5 if args.compact else 2):blocked=True;break
                if not blocked:
                    for reference in reference_poses[qi]:
                        if wire_pair_gap(path,reference)<1:blocked=True;break
                if qi==0 and not blocked and args.space_check:
                    u=(G-A)/np.linalg.norm(G-A);nn=an-u*(an@u);nn/=np.linalg.norm(nn);vv=np.cross(nn,u);smax=np.linalg.norm(G-A)-24
                    belly=md.Manifold.cube([smax,34,14],center=True).transform(np.column_stack([u,vv,nn,A+u*smax/2]))
                    uw=(d.geom_xpos[gid]-d.site_xpos[ids['A']]);uw/=np.linalg.norm(uw);wheel=md.Manifold.cylinder(6,9.5,circular_segments=24,center=True).transform(np.column_stack([uw,np.cross(routed['axis'],uw),routed['axis'],routed['center']]))
                    for bone,part in parts.items():
                        body=part['body'];rb=d.xmat[body].reshape(3,3);rc=d.xmat[cb].reshape(3,3);pb=d.xpos[body]*1000;pcw=d.xpos[cb]*1000
                        shapes=[belly.transform(np.column_stack([rb.T@rc,rb.T@(pcw-pb)])),wheel.transform(np.column_stack([rb.T,-rb.T@pb]))]
                        for shape in shapes:
                            box=np.array(shape.bounding_box())
                            if np.any(box[:3]>=part['hi']) or np.any(box[3:]<=part['lo']):continue
                            if (shape^part['solid']).volume()>.001:blocked=True;break
                        if blocked:break
                if qi>0 and not blocked and args.space_check:
                    uw=(d.geom_xpos[gid]-d.site_xpos[ids['A']]);uw/=np.linalg.norm(uw);wheel=md.Manifold.cylinder(6,9.5,circular_segments=24,center=True).transform(np.column_stack([uw,np.cross(routed['axis'],uw),routed['axis'],routed['center']]))
                    for part in parts.values():
                        rb=d.xmat[part['body']].reshape(3,3);localwheel=wheel.transform(np.column_stack([rb.T,-rb.T@(d.xpos[part['body']]*1000)]));box=np.array(localwheel.bounding_box())
                        if np.any(box[:3]>=part['hi']) or np.any(box[3:]<=part['lo']):continue
                        if (localwheel^part['solid']).volume()>.001:blocked=True;break
                if qi==0:
                    zero=not blocked
                    if blocked:break
                    if not functional:break
                    if old_matrix is not None:
                        replacement=old_matrix.copy()
                        for jj in range(7):
                            plusq=np.zeros(7);plusq[jj]=1e-4;rp2=route_at(plusq);rm2=route_at(-plusq)
                            replacement[index,jj]=-(rp2['D_mm']-rm2['D_mm'])/(2*np.radians(1e-4)) if rp2 and rm2 else np.nan
                        replacement[np.abs(replacement)<1e-4]=0
                        closure=bool(np.all(np.isfinite(replacement)) and np.linalg.matrix_rank(replacement)==7 and linprog(np.zeros(16),A_eq=np.vstack([replacement.T,np.ones(16)]),b_eq=np.r_[np.zeros(7),1],bounds=(1e-6,None),method='highs').success)
                        if not closure:functional=False;break
                if not blocked:clear+=1
            score=(zero and functional,zero,clear,worst,-offset)
            if best is None or score>best['score']:best={'score':score,'A':A.copy(),'G':G.copy(),'B':B.copy(),'afoot':afoot.copy(),'gfoot':gfoot.copy(),'bfoot':bfoot.copy(),'normal':an.copy(),'zero_clear':zero,'coarse_clear':clear,'offset':offset,'L_AG':float(np.linalg.norm(G-A)),'functional':functional,'moment':moment}
            if zero and functional and clear==len(coarse):break
        if best is None:raise RuntimeError('No placement candidates '+rid)
        if args.preserve_closure and (not best['functional'] or not best['zero_clear']):
            for role in ['A','C','B','Aroot','Groot','Broot']:root.find(f'.//site[@name="{rid}_{role}"]').set('pos',cached_root.find(f'.//site[@name="{rid}_{role}"]').get('pos'))
            root.find(f'.//geom[@name="{rid}_guide"]').set('pos',cached_root.find(f'.//geom[@name="{rid}_guide"]').get('pos'));results.append(cached[rid]);m,d=compile_tree(tree);print('No passing replacement; retained',rid,flush=True);continue
        for role,p in [('A',best['A']),('C',best['A']),('B',best['B']),('Aroot',best['afoot']),('Groot',best['gfoot']),('Broot',best['bfoot'])]:root.find(f'.//site[@name="{rid}_{role}"]').set('pos',' '.join(str(x*.001) for x in p))
        root.find(f'.//geom[@name="{rid}_guide"]').set('pos',' '.join(str(x*.001) for x in best['G']))
        record={'id':rid,'carrier':carrier,'insertion_body':target,'belly_fixed_eye_local_mj_mm':best['A'].tolist(),'guide_center_local_mj_mm':best['G'].tolist(),'guide_coordinate_semantics':'fixed incoming tangent/pivot G, not sheave center; R8 sheave center swivels around incoming A-G axis','insertion_eye_local_mj_mm':best['B'].tolist(),'belly_root_local_mj_mm':best['afoot'].tolist(),'guide_root_local_mj_mm':best['gfoot'].tolist(),'insertion_root_local_mj_mm':best['bfoot'].tolist(),'belly_normal_local_mj':best['normal'].tolist(),'available_AG_mm':best['L_AG'],'target_joint_index':target_axis,'target_sign':target_sign,'target_moment_arm_mm':best['moment'],'target_direction_screen_pass':best['functional'],'coarse_steel_zero_clear':best['zero_clear'],'coarse_steel_clear_poses':best['coarse_clear'],'coarse_pose_count':len(coarse),'belly_pose_solved':False,'TPU_volume_verified':False,'hardware_verified':False,'frozen':False}
        results.append(record);print(json.dumps({k:v for k,v in record.items() if not isinstance(v,list)}),flush=True)
        m,d=compile_tree(tree)
        OUT.mkdir(exist_ok=True,parents=True);tree.write(OUT/'psm_module_candidate.xml',encoding='utf-8');(OUT/'psm_module_candidate.json').write_text(json.dumps({'architecture':'A fixed TPU belly end, C moving TPU end and active-wire termination, passive tendon C -> R8 swivelling sheave entry G -> bone eye B; active feed motor -> A -> C. C=A only for corridor search. XML spheres mark G; actual sheave center/arc use swivel_route, not the XML sphere proxy.','modules':results},indent=2))
    tree.write(OUT/'psm_module_candidate.xml',encoding='utf-8');(OUT/'psm_module_candidate.json').write_text(json.dumps({'architecture':'fixed A, moving C, passive tendon via swivelling sheave entry G to B; active motor wire terminates at C; XML G spheres are markers, use analytic swivel_route geometry','modules':results},indent=2))
def oldmujocoid(m,name):return mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,name)
if __name__=='__main__':run()
