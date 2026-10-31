"""Prescribed-path solid probe, including available candidate joint hardware.
Cartilage, muscle deformation, strength and manufacturing remain unvalidated.
"""
from repair_inputs import *
from scipy.spatial.transform import Rotation
import mujoco
DEST=OUT/'joint_model'
WITH_MOUNTS='--muscle-mounts' in sys.argv
WITH_ELBOW='--elbow' in sys.argv or WITH_MOUNTS
WITH_GH='--GH' in sys.argv or WITH_ELBOW
WITH_AC='--AC' in sys.argv or WITH_GH
WITH_SC='--SC' in sys.argv or WITH_AC
SC_DEST=DEST/'SC_connection'
if WITH_SC:DEST=SC_DEST
AC_DEST=SC_DEST/('AC_compact' if '--compact' in sys.argv or WITH_GH else 'AC_connection')
if WITH_AC:DEST=AC_DEST
if WITH_GH:DEST=DEST/('GH_closed' if '--closed' in sys.argv else 'GH_preserve_head' if '--preserve-head' in sys.argv else 'GH_connection')
GH_DEST=DEST
if WITH_ELBOW:DEST=DEST/'elbow_connection'
ELBOW_DEST=DEST
if WITH_MOUNTS:DEST=OUT/'elbow_muscles/mounted_elbow_candidate'

def load(path):
    v,f=read_obj(str(path));s=mf.Manifold(mf.Mesh64((v-SHIFT).astype(np.float64),f.astype(np.uint64)))
    if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid solid '+str(path))
    return s

def controls_to_q(m,c):
    q=m.qpos0.copy()
    for name,i in [('SC',0),('AC',3),('GH',6)]:
        j=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,name+'_spherical');a=int(m.jnt_qposadr[j]);xyzw=Rotation.from_rotvec(np.radians(c[i:i+3])).as_quat();q[a:a+4]=xyzw[[3,0,1,2]]
    q[-1]=np.radians(c[9]);return q

def run():
    # Capture the actual mesh inputs before simulation. XML hashes alone only
    # identify path strings and cannot detect replacement geometry at a path.
    scene=ET.parse(DEST/'arm.xml').getroot();compiler=scene.find('compiler')
    meshdir=compiler.get('meshdir','.') if compiler is not None else '.'
    mesh_evidence={str((DEST/meshdir/mesh.get('file')).resolve()):hashlib.sha256((DEST/meshdir/mesh.get('file')).read_bytes()).hexdigest() for mesh in scene.findall('./asset/mesh')}
    m=mujoco.MjModel.from_xml_path(str(DEST/'arm.xml'));d=mujoco.MjData(m);definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'));inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']};owners=definition['bone_owner'];centres={k:np.array(v)-SHIFT for k,v in definition['centres_source_mm'].items()};moving={}
    for fid in owners:
        v,f=read_obj(str(ROOT/byid[fid]['source']));v,f=weld_stl(v[f]);s=mf.Manifold(mf.Mesh64((v-SHIFT).astype(np.float64),f.astype(np.uint64)))
        if s.status()!=mf.Error.NoError:raise RuntimeError(fid+' not closed')
        moving[fid]=s
    fixed={r['id']:load(OUT/'repaired_supports'/(r['id']+'.obj')) for r in support['models'] if r['kind']=='bone_reference'}
    if WITH_SC:
        assembly=json.loads((SC_DEST/'assembly.json').read_text(encoding='utf8'))
        for name in assembly['components']:
            dest=moving if name=='clavicle_SC_ball' else fixed
            v,f=read_obj(str(SC_DEST/(name+'.obj')));s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid SC part '+name)
            dest[{'clavicle_SC_ball':'FJ3362','manubrium_socket_base':'FJ3290'}.get(name,name)]=s
    if WITH_AC:
        assembly=json.loads((AC_DEST/'assembly.json').read_text(encoding='utf8'))
        for name,owner in assembly['part_owners'].items():
            fid={'clavicle_SC_AC_balls':'FJ3362','scapula_AC_socket':'FJ3384'}.get(name,name)
            v,f=read_obj(str(AC_DEST/(name+'.obj')))
            s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid AC part '+name)
            moving[fid]=s;owners[fid]=owner
    if WITH_GH:
        assembly=json.loads((GH_DEST/'assembly.json').read_text(encoding='utf8'))
        for name,owner in assembly['part_owners'].items():
            fid=assembly['bone_aliases'].get(name,name);v,f=read_obj(str(GH_DEST/(name+'.obj')))
            s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid GH part '+name)
            moving[fid]=s;owners[fid]=owner
    if WITH_ELBOW:
        assembly=json.loads((ELBOW_DEST/'assembly.json').read_text(encoding='utf8'))
        for fid in assembly['removed_geoms']:moving.pop(fid);owners.pop(fid)
        for name,owner in assembly['part_owners'].items():
            fid=assembly['bone_aliases'].get(name,name);v,f=read_obj(str(ELBOW_DEST/(name+'.obj')))
            s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid elbow part '+name)
            moving[fid]=s;owners[fid]=owner
    if WITH_MOUNTS:
        assembly=json.loads((DEST/'assembly.json').read_text(encoding='utf8'))
        for r in assembly['instances']:
            if r['owner'] in ['soft_unvalidated','distal_unvalidated']:continue
            v,f=read_obj(str(ROOT/r['file']));s=mf.Manifold(mf.Mesh64((v+np.array(r['translation_display_mm'])).astype(np.float64),f.astype(np.uint64)))
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid M2 instance '+r['key'])
            if r['owner']=='world':fixed[r['key']]=s
            else:moving[r['key']]=s;owners[r['key']]=r['owner']
    fboxes={k:np.array(s.bounding_box()) for k,s in fixed.items()};pairs=list(itertools.combinations(moving,2));cache={};bids={n:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,n) for n in centres};GH=centres['GH'];E=centres['elbow']
    rigid_pair_cache={}
    def pose(c,refresh=False):
        key=tuple(np.round(c,8))
        if key in cache and not refresh:return cache[key]
        d.qpos[:]=controls_to_q(m,c);mujoco.mj_forward(m,d);world={}
        for fid,s in moving.items():
            owner=owners[fid];b=bids[owner];R=d.xmat[b].reshape(3,3);world[fid]=s.transform(np.c_[R,d.xpos[b]*1000-R@centres[owner]])
        hits=[];gap=3.;closest=[];boxes={k:np.array(s.bounding_box()) for k,s in world.items()}
        def check(a,sa,ba,b,sb,bb,rigid=False):
            nonlocal gap,closest
            if rigid and (a,b) in rigid_pair_cache:
                vol,g=rigid_pair_cache[a,b]
                if vol>.001:hits.append({'pair':[a,b],'volume_mm3':vol})
                if g<gap:gap=g;closest=[a,b]
                return
            delta=np.maximum(ba[:3]-bb[3:],bb[:3]-ba[3:])
            if np.linalg.norm(np.maximum(delta,0))>=3:return
            vol=max(0.,float((sa^sb).volume())) if np.all(delta<0) else 0.
            if vol>.001:hits.append({'pair':[a,b],'volume_mm3':vol})
            g=float(sa.min_gap(sb,3))
            if rigid:rigid_pair_cache[a,b]=(vol,g)
            if g<gap:gap=g;closest=[a,b]
        for a,b in pairs:check(a,world[a],boxes[a],b,world[b],boxes[b],owners[a]==owners[b])
        for a,sa in world.items():
            for b,sb in fixed.items():check(a,sa,boxes[a],b,sb,fboxes[b])
        R=d.xmat[bids['GH']].reshape(3,3);elevation=float(np.degrees(np.arccos(np.clip((R@(E-GH))@(E-GH)/np.linalg.norm(E-GH)**2,-1,1))))
        row={'controls_deg':list(c),'qpos':d.qpos.tolist(),'hits':hits,'min_bone_gap_mm_capped_3':gap,'closest_pair':closest,'humerus_elevation_from_zero_deg':elevation};cache[key]=row;return row
    def c(**kwargs):
        x=np.zeros(10)
        for key,value in kwargs.items():x[int(key[1:])]=value
        return x
    long=(E-GH)/np.linalg.norm(E-GH);flex_axis=np.cross(long,[0,-1,0]);flex_axis/=np.linalg.norm(flex_axis);abd_axis=np.cross(long,[-1,0,0]);abd_axis/=np.linalg.norm(abd_axis)
    coupled=c(a1=15,a4=30);coupled[6:9]=abd_axis*120
    tasks=[('elbow_flexion',c(a9=150)),('GH_flexion',np.r_[np.zeros(6),flex_axis*160,0]),('GH_abduction',np.r_[np.zeros(6),abd_axis*160,0]),('GH_axial_positive',np.r_[np.zeros(6),long*90,0]),('GH_axial_negative',np.r_[np.zeros(6),-long*90,0]),('SC_elevation',c(a1=30)),('SC_protraction',c(a2=30)),('SC_retraction',c(a2=-30)),('SC_rotation_X',c(a0=30)),('AC_upward',c(a4=50)),('AC_rotation_X',c(a3=30)),('AC_rotation_Z',c(a5=30)),('coordinated_raise',coupled)]
    if '--only-probes' in sys.argv:
        requested=set(sys.argv[sys.argv.index('--only-probes')+1].split(','));tasks=[(name,end) for name,end in tasks if name in requested]
        if not tasks:raise RuntimeError('No matching probe names')
    def angular_distance(a,b):return max(np.linalg.norm(a[i:i+3]-b[i:i+3]) for i in [0,3,6]) if abs(a[9]-b[9])==0 else max(abs(a[9]-b[9]),*(np.linalg.norm(a[i:i+3]-b[i:i+3]) for i in [0,3,6]))
    rows=[];paths={};neutral=pose(np.zeros(10))
    neutral_report={'scope':'all included solid pairs at zero pose; not a motion pass','moving_count':len(moving),'fixed_count':len(fixed),'neutral':neutral}
    (DEST/'neutral_check.json').write_text(json.dumps(neutral_report,ensure_ascii=False,indent=2),encoding='utf8')
    print('neutral solid hits',neutral['hits'],flush=True)
    if '--neutral-only' in sys.argv or neutral['hits']:return
    for name,end in tasks:
        steps=max(1,int(np.ceil(angular_distance(end,np.zeros(10))-1e-9)));forward=[];first=None
        for i in range(steps+1):
            p=pose(end*i/steps);forward.append(p)
            if first is None and p['hits']:first=i
        safe=1. if first is None else 0.;boundary=None
        if first is not None:
            boundary=forward[first]
            if first>0:
                lo=(first-1)/steps;hi=first/steps
                while angular_distance(end*hi,end*lo)>.1:
                    mid=(lo+hi)/2
                    if pose(end*mid)['hits']:hi=mid
                    else:lo=mid
                safe=lo;boundary=pose(end*hi)
        replay=[pose(np.array(p['controls_deg']),True) for p in reversed(forward)];accepted=pose(end*safe)
        row={'probe':name,'requested_controls_deg':end.tolist(),'continuous_prefix_fraction':safe,'last_clear_controls_deg':accepted['controls_deg'],'actual_humerus_elevation_deg':accepted['humerus_elevation_from_zero_deg'],'first_limit':boundary['hits'] if boundary else [],'forward_samples':len(forward),'return_samples':len(replay),'return_colliding_samples':sum(bool(p['hits']) for p in replay),'full_prescribed_path_bone_clear':first is None,'engineering_joint_or_drive_pass':False};rows.append(row);paths[name]=forward;print(name,'prefix',round(safe,4),'limit',row['first_limit'][:3],flush=True)
    report={'scope':'current independent 10DOF candidate, five moving bones vs fixed bony references and included SC hardware; no adjacent moving bone pairs skipped','fixed_reference_count':len(fixed),'moving_bone_count':len(moving),'angular_step_max_deg':1,'first_boundary_deg':.1,'volume_threshold_mm3':.001,'gap_cap_mm':3,'neutral':neutral,'actions':rows,'cartilage_policy':'7 cartilage references excluded from rigid bone acceptance; compliance/contacts not validated','not_included':['AC/GH/elbow retention hardware','muscle deformation','actuation','manufacturing tolerance and loads'],'source_hashes':{r['id']:r['sha256'] for r in inv['models']},'frozen':False};(DEST/'bone_motion_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');(DEST/'motion_paths.json').write_text(json.dumps(paths,indent=2),encoding='utf8')

    report.update(scope='Independent 10DOF prescribed-path solid audit; all moving pairs included. Same-owner rigid pair results reused by rigid-transform invariance.',moving_bone_count=5,moving_component_count=len(moving),moving_owners=owners,included_physical_connections=['SC','AC'] if WITH_AC else ['SC'] if WITH_SC else [],not_included=['GH/elbow retention hardware','muscle deformation','actuation','manufacturing tolerance and loads'] if WITH_AC else ['AC/GH/elbow retention hardware','muscle deformation','actuation','manufacturing tolerance and loads'])
    if WITH_GH:report.update(included_physical_connections=['SC','AC','GH'],not_included=['elbow retention hardware','muscle deformation','actuation','manufacturing tolerance and loads'])
    if WITH_ELBOW:report.update(moving_bone_count=4,source_bone_count=5,included_physical_connections=['SC','AC','GH','elbow'],not_included=['muscle deformation','actuation','manufacturing tolerance and loads'],forearm_policy='Source ulna and radius joined by actual material bridges in FJ3391; source radius material retained, no pronation/supination')
    if WITH_MOUNTS:report.update(scope='M2 known rigid proximal grounding only; all included moving pairs checked, same-owner invariance reused. Does not establish whole-module motion.',not_included=['TPU bodies and free distal ends (22 components), no compliant route yet','passive tendon/rotor transmission','actuation, manufacturing tolerance and loads'],rigid_motion_scope=json.loads((DEST/'rigid_motion_scope.json').read_text(encoding='utf8')))
    report['complete_prescribed_suite']=len(tasks)==13
    report['input_evidence_sha256']={name:hashlib.sha256((DEST/name).read_bytes()).hexdigest() for name in ['arm.xml','assembly.json'] if (DEST/name).exists()}
    for path,digest in mesh_evidence.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise RuntimeError('Geometry changed during audit '+path)
    report['mesh_evidence_sha256']=mesh_evidence
    (DEST/'bone_motion_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':run()
