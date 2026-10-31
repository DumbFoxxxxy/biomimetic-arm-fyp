"""Full-arm exact STL audit after E1 integration, including all E1 accessories.
Explicitly excludes uninstalled SC/AC/GH reference primitives from acceptance.
Reports that exclusion as unresolved hardware, never whole-arm PASS.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build import *
import mujoco

DEST=OUT/'elbow_connection'
if '--mounts' in sys.argv:DEST=DEST/'muscle_mounts'
if '--fittings' in sys.argv:DEST=DEST/'muscle_mounts/end_fittings'

def run():
    m=_load_project_model(str(DEST/'arm.xml'));d=mujoco.MjData(m)
    frames=json.loads((OUT/'body_frames.json').read_text(encoding='utf8'))['origins_mm'];lm=json.loads((OUT/'source_audit.json').read_text(encoding='utf8'))['landmarks_mm'];E=np.array(lm['E']);GH=np.array(lm['GH'])
    tree=ET.parse(DEST/'arm.xml');root=tree.getroot();source={x.get('name'):x.get('file') for x in root.findall('./asset/mesh')}
    names=['clavicle','scapula','humerus','ulna','elbow_pin','cotter_pin','thrust_washer_-1','thrust_washer_1'];parts={};owners={}
    if '--fittings' in sys.argv:
        for prefix in ['elbow_flexor','elbow_extensor']:
            names.extend(prefix+'_'+suffix for suffix in ['A_ferrule','A_spacer_-1','A_spacer_1'])
    for n in names:
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,n);bid=int(m.geom_bodyid[gid]);owner=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,bid);geom=root.find(f'.//geom[@name="{n}"]');p=(DEST/source[geom.get('mesh')]).resolve();v,f=read_stl(str(p));s=mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))
        if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid mesh '+n)
        parts[n]=s;owners[n]=(bid,np.array(frames[owner]))
    pairs=list(itertools.combinations(names,2));cache={}
    j={mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,i):int(m.jnt_qposadr[i]) for i in range(m.njnt)}
    def pose(q,refresh=False):
        key=tuple(np.round(q,8))
        if key in cache and not refresh:return cache[key]
        d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);world={}
        for n,s in parts.items():
            bid,o=owners[n];R=d.xmat[bid].reshape(3,3);world[n]=s.transform(np.c_[R,d.xpos[bid]*1000-R@o])
        hits=[];gap=5.;closest=[]
        for an,bn in pairs:
            aa=np.array(world[an].bounding_box());bb=np.array(world[bn].bounding_box());delta=np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:]);lower=np.linalg.norm(np.maximum(delta,0))
            if lower>=5:continue
            volume=float((world[an]^world[bn]).volume()) if np.all(delta<0) else 0.
            if volume>.001:hits.append({'pair':[an,bn],'volume_mm3':volume})
            # Same-body clearances do not limit motion but still checked for penetration.
            if owners[an][0]!=owners[bn][0]:
                g=float(world[an].min_gap(world[bn],5))
                if g<gap:gap=g;closest=[an,bn]
        bid=owners['humerus'][0];vec=d.xmat[bid].reshape(3,3)@(E-GH);elev=np.degrees(np.arccos(np.clip(unit(vec)@unit(E-GH),-1,1)))
        row={'q_deg':list(q),'hits':hits,'minimum_component_gap_mm_capped_5':gap,'closest_pair':closest,'actual_humerus_elevation_deg':float(elev)};cache[key]=row;return row
    def q(**values):
        x=np.zeros(m.nq)
        for n,v in values.items():x[j[n]]=v
        return x
    zero=q();tasks=[('K1',[zero,q(elbow=145)]),('K2',[zero,q(gh_flex=120)]),('K3',[zero,q(gh_abd=120)]),('K4',[zero,q(gh_rot=70)]),('K5',[zero,q(gh_rot=-90)]),('K6',[zero,q(ac_up=50)]),('K7',[zero,q(ac_up=50,gh_abd=120)]),('RAISE_THEN_FLEX_ELBOW',[zero,q(gh_flex=96),q(gh_flex=96,elbow=145),q(elbow=145),zero]),('ABDUCT_THEN_FLEX_ELBOW',[zero,q(gh_abd=96),q(gh_abd=96,elbow=145),q(elbow=145),zero]),('HUG',[zero,q(gh_flex=25,gh_abd=8,gh_rot=-20),q(gh_flex=25,gh_abd=8,gh_rot=-20,elbow=116)]),('SC_FORWARD',[zero,q(sc_protract=5)]),('SC_BACK',[zero,q(sc_protract=-5)])]
    rows=[];paths={}
    for name,knots in tasks:
        samples=[];first=None;safe=zero.copy();first_index=None
        for seg,(a,b) in enumerate(zip(knots[:-1],knots[1:])):
            steps=max(1,int(np.ceil(max(abs(b-a)))))
            for i in range(steps+1):
                if seg and i==0:continue
                value=a+(b-a)*i/steps;p=pose(value);samples.append(p)
                if first is None:
                    if p['hits']:
                        first=p;first_index=len(samples)-1;lo=safe.copy();hi=value.copy()
                        while max(abs(hi-lo))>.1:
                            mid=(hi+lo)/2
                            if pose(mid)['hits']:hi=mid
                            else:lo=mid
                        safe=lo;first=pose(hi)
                    else:safe=value.copy()
        reverse=[pose(np.array(p['q_deg']),True) for p in reversed(samples)]
        rows.append({'action':name,'continuous_prefix_last_clear':pose(safe),'first_limit':first,'first_collision_sample_index':first_index,'forward_samples':len(samples),'return_samples':len(reverse),'return_colliding_samples':sum(bool(p['hits']) for p in reverse),'whole_path_component_clear':first is None,'whole_arm_accepted':False})
        paths[name]=samples;print(name,'prefix',np.round(safe,2).tolist(),'limit',first['hits'] if first else [],flush=True)
    neutral_muscles=[]
    for n in ['elbow_flexor','elbow_extensor']:
        v,f=read_stl(str(OUT/'assets'/(n+'.stl')));s=mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)));hits=[]
        for pn,p in parts.items():
            vol=float((s^p).volume())
            if vol>.001:hits.append({'part':pn,'volume_mm3':vol})
        neutral_muscles.append({'muscle':n,'solid_neutral_overlap':hits,'compliant_acceptance':False})
    dump(DEST/'full_arm_audit.json',{'status':'E1 integrated geometry audit, not complete engineering arm acceptance','included_parts':names,'all_pair_count':len(pairs),'angular_step_max_deg':1,'first_boundary_resolution_deg':.1,'volume_threshold_mm3':.001,'gap_cap_mm':5,'unresolved_excluded_reference_geometry':['SC_pin','AC_pin','AC_bearing_reference','GH_ball_reference','frame_mount'],'qualification':'SC/AC/GH supports and real frame still uninstalled. They cannot receive acceptance from this audit. All included part pairs checked including same-owner pairs.','actions':rows,'neutral_muscles':neutral_muscles});dump(DEST/'full_arm_motion_paths.json',paths)

if __name__=='__main__':run()
