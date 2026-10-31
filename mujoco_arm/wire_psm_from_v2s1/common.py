
# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys,json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'mujoco_arm/skeleton_from_v2'),str(ROOT/'mujoco_arm/blade_motion_design')]
from audit_baseline import local_solid
from geometry_distance import segment_mesh_distance
BONES=['thorax_geom','clavicle_visual','scapula_mesh','humerus_shaft','ulna_shaft']
OUT=ROOT/'mujoco_arm/exports/wire_psm_from_v2s1'
BASE=ROOT/'mujoco_arm/exports/skeleton_from_v2/V2-S1'

def template():
    tree=ET.parse(ROOT/'mujoco_arm/cable_design_v2/arm_cable_design_v2.xml');r=tree.getroot();r.set('model','V2S1_wire_PSM_interface_candidate')
    for tag in ['actuator','sensor','contact','equality']:
        x=r.find(tag)
        if x is not None:r.remove(x)
    r.find('option').set('gravity','0 0 0')
    world=r.find('worldbody');thorax=r.find('.//body[@name="thorax"]')
    # Old S05/S06 fixed sites were children of world, although the brief
    # requires attachments on the arm/thorax. Same zero frame, physical owner.
    routed={n.get('site') for t in r.find('tendon') for n in t if n.tag=='site'}
    for site in list(world.findall('site')):
        if site.get('name') in routed:world.remove(site);thorax.append(site)
    for parent in r.iter():
        for child in list(parent):
            if child.tag=='geom' and child.get('name') is not None and child.get('name') not in BONES and not child.get('name','').startswith('V2_guide'):parent.remove(child)
    for n in BONES:
        g=r.find(f'.//geom[@name="{n}"]');asset=r.find(f'.//asset/mesh[@name="{g.get("mesh")}"]')
        asset.set('file',str(BASE/'assets'/(n+'.stl')));g.set('pos','0 0 0')
    used={g.get('mesh') for g in r.findall('.//geom') if g.get('mesh')}
    for x in list(r.find('asset')):
        if x.tag=='mesh' and x.get('name') not in used:r.find('asset').remove(x)
    return tree

def compile_tree(tree):
    compiler=tree.getroot().find('compiler')
    if compiler is not None and compiler.get('meshdir')=='.':compiler.set('meshdir',str(BASE))
    m=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'))
    # Positions are optimized after compilation. Identity-frame fast paths
    # would otherwise ignore position changes for sites/guides born at zero.
    m.site_sameframe[:]=0;m.geom_sameframe[:]=0
    d=mujoco.MjData(m);mujoco.mj_forward(m,d);return m,d

def solids(m,d):
    out={}
    for n in BONES:
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,n);s,_=local_solid(m,d,g);mesh=s.to_mesh();tri=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
        out[n]={'solid':s,'body':int(m.geom_bodyid[g]),'tri':tri,'lo':tri.min((0,1)),'hi':tri.max((0,1)),'tlo':tri.min(1),'thi':tri.max(1)}
    return out

def project(p,tri):
    a,b,c=tri[:,0],tri[:,1],tri[:,2];normal=np.cross(b-a,c-a);normal/=np.maximum(np.linalg.norm(normal,axis=1),1e-20)[:,None]
    pr=p-np.sum((p-a)*normal,axis=1)[:,None]*normal
    v0=b-a;v1=c-a;v2=pr-a;d00=np.sum(v0*v0,axis=1);d01=np.sum(v0*v1,axis=1);d11=np.sum(v1*v1,axis=1);d20=np.sum(v2*v0,axis=1);d21=np.sum(v2*v1,axis=1);den=d00*d11-d01*d01
    u=(d11*d20-d01*d21)/np.maximum(den,1e-20);v=(d00*d21-d01*d20)/np.maximum(den,1e-20);valid=(u>=0)&(v>=0)&(u+v<=1)&(den>1e-12)
    candidates=[pr];errors=[np.where(valid,np.sum((pr-p)**2,axis=1),np.inf)]
    for x,y in [(a,b),(b,c),(c,a)]:
        delta=y-x;t=np.clip(np.sum((p-x)*delta,axis=1)/np.maximum(np.sum(delta*delta,axis=1),1e-20),0,1);z=x+t[:,None]*delta;candidates.append(z);errors.append(np.sum((z-p)**2,axis=1))
    flat=np.array(errors);i,j=np.unravel_index(np.argmin(flat),flat.shape);return candidates[i][j],normal[j],float(np.sqrt(flat[i,j]))

def paths(m,d):
    xyz=d.wrap_xpos.reshape(-1,3)*1000;obj=d.wrap_obj.ravel();result=[]
    for i in range(m.ntendon):
        a=int(d.ten_wrapadr[i]);n=int(d.ten_wrapnum[i]);p=xyz[a:a+n];o=obj[a:a+n]
        if n in [2,3]:result.append(p.copy());continue
        if n!=4 or o[1]!=o[2]:raise RuntimeError('Unsupported wrap topology: '+str((i,n,o)))
        p0,p1,p2,p3=p;g=int(o[1]);normal=d.geom_xmat[g].reshape(3,3)[:,2];origin=d.geom_xpos[g]*1000
        if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_SPHERE:
            r1=p1-origin;r2=p2-origin;normal=np.cross(r1,r2);normal/=max(np.linalg.norm(normal),1e-20)
            theta=np.arccos(np.clip(r1@r2/max(np.linalg.norm(r1)*np.linalg.norm(r2),1e-20),-1,1));t=np.linspace(0,theta,max(2,int(np.ceil(theta/np.radians(2)))+1))
            arc=origin+r1[None,:]*np.cos(t)[:,None]+np.cross(normal,r1)[None,:]*np.sin(t)[:,None]
            result.append(np.vstack([p0,arc,p3]));continue
        z1=normal@(p1-origin);z2=normal@(p2-origin);r1=p1-origin-normal*z1;r2=p2-origin-normal*z2
        angle=np.arctan2(normal@np.cross(r1,r2),r1@r2);radius=m.geom_size[g,0]*1000
        length=d.ten_length[i]*1000-np.linalg.norm(p1-p0)-np.linalg.norm(p3-p2)
        theta=min([angle,angle+2*np.pi,angle-2*np.pi],key=lambda x:abs(np.hypot(x*radius,z2-z1)-length))
        t=np.linspace(0,1,max(2,int(np.ceil(abs(theta)/np.radians(2)))+1));rr=r1[None,:]*np.cos(t*theta)[:,None]+np.cross(normal,r1)[None,:]*np.sin(t*theta)[:,None]
        result.append(np.vstack([p0,origin+normal*(z1+(z2-z1)*t)[:,None]+rr,p3]))
    return result

def local_gap(path,part,margin=2.5):
    lo=part['lo'];hi=part['hi'];minimum=1e6
    if np.any(path.min(0)-margin>hi) or np.any(path.max(0)+margin<lo):return minimum
    for a,b in zip(path[:-1],path[1:]):
        selected=np.all(part['tlo']<=np.maximum(a,b)+margin,axis=1)&np.all(part['thi']>=np.minimum(a,b)-margin,axis=1)
        if selected.any():minimum=min(minimum,segment_mesh_distance(a,b,part['tri'][selected])-.5)
        mid=(a+b)/2
        if np.all(mid>lo) and np.all(mid<hi) and len(part['solid'].ray_cast((lo-[123,151,179]).tolist(),mid.tolist()))%2:minimum=min(minimum,-.5)
        if minimum<0:break
    return minimum

def pose_library():return json.loads((BASE/'motion_paths.json').read_text())
def design(x):return np.asarray(x)[...,[0,2,1]]*np.array([1,1,-1])

def swivel_route(A,G,B,radius=8.):
    """Swivelling sheave, fixed incoming tangent G and A-G entry axis.

    The wheel center is offset R from G, rotates about the incoming axis,
    and is not a fictitious fixed spherical obstacle. Return exact arc length.
    """
    u=np.asarray(G)-A;u/=max(np.linalg.norm(u),1e-20);b=np.asarray(B)-G;bu=float(b@u);transverse=b-bu*u;length=np.linalg.norm(transverse)
    if length<1e-8:
        helper=np.array([1.,0,0]) if abs(u[0])<.8 else np.array([0.,1.,0]);v=np.cross(u,helper);v/=np.linalg.norm(v)
    else:v=transverse/length
    bv=float(b@v)-radius;d2=bu*bu+bv*bv
    if d2<=radius*radius+1e-8:return None
    base=radius*radius/d2*np.array([bu,bv]);perp=np.array([-bv,bu]);side=radius*np.sqrt(d2-radius*radius)/d2*perp;choices=[]
    for point in [base+side,base-side]:
        theta=float(np.arctan2(point[0],-point[1])%(2*np.pi))
        if abs(theta-2*np.pi)<1e-8:theta=0.
        tangent=np.array([np.cos(theta),np.sin(theta)]);out=np.array([bu,bv])-point
        if out@tangent>0:choices.append(theta)
    if not choices:return None
    theta=min(choices);O=G+radius*v;n=np.cross(u,v);t=np.linspace(0,theta,max(2,int(np.ceil(theta/np.radians(2)))+1));arc=O+radius*np.sin(t)[:,None]*u-radius*np.cos(t)[:,None]*v
    return {'path':np.vstack([A,arc,B]),'D_mm':float(radius*theta+np.sqrt(d2-radius*radius)),'center':O,'axis':n,'wrap_deg':float(np.degrees(theta))}
