"""Actual material edits for the open concentric prototype, not collision masks."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import xml.etree.ElementTree as ET
import json
import numpy as np
import mujoco
from shoulder_model import v8,EXPORT,inputs,SOURCE,add_connector
from search_gh_anchors import A,B,C
from stl_io import write_binary_stl
import sys
sys.path.insert(0,str(SOURCE/'vendor'))
import manifold3d as md

def cyl(axis,center,length,radius,bore=0.):
    helper=C if abs(np.dot(axis,C))<.9 else B
    u=np.cross(helper,axis);u/=np.linalg.norm(u);vv=np.cross(axis,u)
    solid=md.Manifold.cylinder(length,radius,circular_segments=48,center=True)
    if bore:solid=solid-md.Manifold.cylinder(length+4,bore,circular_segments=48,center=True)
    return solid.transform(np.column_stack([np.column_stack([u,vv,axis]),center]))

def rod(a,b,radius):
    a=np.array(a);b=np.array(b);delta=b-a
    return cyl(delta/np.linalg.norm(delta),(a+b)/2,np.linalg.norm(delta),radius)+md.Manifold.sphere(radius,circular_segments=24).translate(a)+md.Manifold.sphere(radius,circular_segments=24).translate(b)

def mesh_solid(model,data,name):
    from audit_geometry import mesh_local
    gid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,name)
    _,tri=mesh_local(model,data,gid);tri=tri*1000@v8.DESIGN_TO_MJ
    verts,faces=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True)
    return md.Manifold(md.Mesh(np.ascontiguousarray(verts,dtype=np.float32),np.ascontiguousarray(faces.reshape(-1,3),dtype=np.uint32)))

def install(root,owner,name,solid,color='.7 .72 .77 1'):
    mesh=solid.to_mesh();p=np.asarray(mesh.vert_properties)[:,:3]
    path=EXPORT/f'{name}_material_v2.stl'
    write_binary_stl(path,(p@v8.DESIGN_TO_MJ.T*.001)[np.asarray(mesh.tri_verts)],'Actual machined prototype geometry; body-local MuJoCo metres')
    asset=root.find('asset');old=asset.find(f'mesh[@name="{name}_material"]')
    if old is not None:asset.remove(old)
    ET.SubElement(asset,'mesh',{'name':name+'_material','file':str(path)})
    old=owner.find(f'geom[@name="{name}"]')
    if old is not None:owner.remove(old)
    ET.SubElement(owner,'geom',{'name':name,'type':'mesh','mesh':name+'_material','pos':'0 0 0','rgba':color,'contype':'0','conaffinity':'0','density':'0'})

def refine(root):
    cfg,fit,anchors=inputs();ac=np.array(cfg['geometry_mm']['AC']);sc=np.array(cfg['geometry_mm']['SC'])
    from gimbal_arm_model import GH_SHIFT
    gh=np.array([159.06,187.54,7.30])+GH_SHIFT
    model=_load_project_xml(ET.tostring(root,encoding='unicode'));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    scap=root.find('.//body[@name="scapula"]');clav=root.find('.//body[@name="clavicle_body"]')
    # Sweep clearance for the medial GF rotor, including a 2 mm machining gap.
    void=cyl(B,gh-ac-50*B,40.,20.)
    before=mesh_solid(model,data,'scapula_mesh');after=before-void
    # The fixed GF inlet-to-drum tangent is fixed in the scapula frame.
    # Cut its actual free-running corridor instead of permitting penetration.
    tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'GF_cable_E02');adr=int(data.ten_wrapadr[tid]);world=data.wrap_xpos.reshape(-1,3)[adr:adr+2]
    sid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'scapula');local=(world-data.xpos[sid])@data.xmat[sid].reshape(3,3)*1000@v8.DESIGN_TO_MJ
    after=after-rod(local[0],local[1],1.2)
    install(root,scap,'scapula_mesh',after)
    # GR uses fixed Bowden ferrules, not filled unused guide cylinders.
    # Each ferrule's bore is aligned with its actual fixed tangent segment.
    carrier=root.find('.//body[@name="gh_abd"]');bid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'gh_abd')
    for node in list(carrier):
        name=node.get('name','')
        if node.tag=='geom' and (name.startswith('GR_entry_') or name.startswith('GR_E01_guide_') or name.startswith('GR_E02_guide_')):carrier.remove(node)
    for rid in ['E01','E02']:
        tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'GR_cable_'+rid);adr=int(data.ten_wrapadr[tid]);world=data.wrap_xpos.reshape(-1,3)[adr:adr+2]
        p=(world-data.xpos[bid])@data.xmat[bid].reshape(3,3)*1000@v8.DESIGN_TO_MJ;direction=p[1]-p[0];direction/=np.linalg.norm(direction);centre=p[0]-6*direction
        install(root,carrier,'GR_ferrule_'+rid,cyl(direction,centre,8.,2.8,.8),'.85 .35 .2 1')
        radial=centre-A*np.dot(centre,A);radial/=np.linalg.norm(radial)
        root_point=-10*A+25*radial;via=-10*A+np.linalg.norm(centre-A*np.dot(centre,A))*radial;attach=centre-3*A
        for i,(a,b) in enumerate([(root_point,via),(via,attach)]):add_connector(carrier,f'GR_ferrule_support_{rid}_{i}',a,b,1.8)
    # A complete ball/socket interface replaces mutually intersecting webs.
    # The cup opens anteriorly to admit the clavicular stem.
    before_clav=mesh_solid(model,data,'clavicle_mesh')
    after_clav=before_clav-md.Manifold.sphere(9.,circular_segments=48).translate(ac-sc)-md.Manifold.sphere(24.,circular_segments=48)
    install(root,clav,'clavicle_mesh',after_clav)
    for owner,name in [(clav,'ac_clavicle_web'),(scap,'ac_scapula_web')]:
        old=owner.find(f'geom[@name="{name}"]')
        if old is not None:owner.remove(old)
    ball=md.Manifold.sphere(4.,circular_segments=48).translate(ac-sc)
    # Fuse the ball with its stem and remove the cup-facing obsolete neck.
    install(root,clav,'ac_ball',ball)
    add_connector(clav,'ac_ball_stem',ac-sc,ac-sc-8*C,2.)
    vertices=np.asarray(after_clav.to_mesh().vert_properties)[:,:3];aim=(ac-sc)/np.linalg.norm(ac-sc)*27.-8*C
    attach=vertices[np.argmin(np.linalg.norm(vertices-aim,axis=1))]
    add_connector(clav,'sc_root_bridge',np.zeros(3),attach,2.8)
    cup=md.Manifold.sphere(6.5,circular_segments=48)-md.Manifold.sphere(4.6,circular_segments=48)
    cup=cup-md.Manifold.cube((40.,40.,20.),center=True).translate([0.,0.,-10.])
    install(root,scap,'ac_socket',cup)
    vertices=np.asarray(after.to_mesh().vert_properties)[:,:3];aim=12*v8.source_basis(cfg)[2]
    attach=vertices[np.argmin(np.linalg.norm(vertices-aim,axis=1))]
    add_connector(scap,'ac_socket_web',6*v8.source_basis(cfg)[2],attach,2.)
    # Carry S06 round the machined GF pocket; retain its exact cable endpoint.
    row=next(r for r in anchors['routes'] if r['id']=='S06');base=np.array([159.06,187.54,7.30])-ac+12*v8.source_basis(cfg)[2];eye=np.array(row['eye_local_ac_mm'])
    old=scap.find('geom[@name="eye_web_S06"]');scap.remove(old)
    via1=base+30*C;via2=eye+25*C
    for i,(a,b) in enumerate([(base,via1),(via1,via2),(via2,eye)]):add_connector(scap,f'eye_web_S06_{i}',a,b,1.5)
    routes_path=EXPORT/'full_arm_shoulder_routes_candidate.json'
    if routes_path.exists():
        extra_routes=json.loads(routes_path.read_text())
        if extra_routes['pass']:
            for r in extra_routes['routes']:
                if 'additional_bracket_offset_local_mm' in r:
                    end=np.array(r['eye_local_ac_mm']);start=end-np.array(r['additional_bracket_offset_local_mm'])
                    if r['id']=='S09':
                        points=[start,start+20*C,end+20*C,end]
                        planned=EXPORT/'s09_bracket_path.json'
                        if planned.exists():
                            p=json.loads(planned.read_text())
                            if p['pass'] and np.allclose(p['path_local_ac_mm'][0],start) and np.allclose(p['path_local_ac_mm'][-1],end):points=[np.array(x) for x in p['path_local_ac_mm']]
                        for i,(a,b) in enumerate(zip(points[:-1],points[1:])):add_connector(scap,f'eye_web_S09_extension_{i}',a,b,1.5)
                    else:add_connector(scap,f'eye_web_{r["id"]}_extension',start,end,1.5)
    # Replace the overlapping legacy elbow ends with axially separated links.
    # The joint centre and all commanded angles remain unchanged.
    hum=root.find('.//body[@name="humerus"]');ulna=root.find('.//body[@name="ulna"]')
    elbow=np.array([184.34,-101.35,7.30])-np.array([159.06,187.54,7.30])
    upper=mesh_solid(model,data,'humerus_shaft')-md.Manifold.sphere(65.,circular_segments=48).translate(elbow)
    fore=mesh_solid(model,data,'ulna_shaft')-md.Manifold.sphere(65.,circular_segments=48)
    parent_end=rod(elbow-70*A,elbow-50*A-12*B,6.)+rod(elbow-50*A-12*B,elbow-8*A-12*B,6.)+cyl(B,elbow-12*B,8.,10.)
    parent_end=parent_end-cyl(B,elbow,100.,3.5)
    child_end=rod(70*A,50*A+12*B,6.)+rod(50*A+12*B,8*A+12*B,6.)+cyl(B,12*B,8.,10.)+cyl(B,7*B,50.,3.)
    install(root,hum,'humerus_shaft',upper+parent_end,'.82 .8 .72 1')
    install(root,ulna,'ulna_shaft',fore+child_end,'.82 .8 .72 1')
    old=hum.find('geom[@name="E_holder_root"]')
    if old is not None:hum.remove(old)
    add_connector(hum,'E_holder_root',elbow-70*A,elbow-35*A+20*B+25*C,2.)
    (EXPORT/'gimbal_material_changes.json').write_text(json.dumps({'scope':'Prototype material edits; original CAD and shoulder-only checkpoint retained',
        'scapula_removed_mm3':before.volume()-after.volume(),'scapula_components':len(after.decompose()),
        'clavicle_removed_mm3':before_clav.volume()-after_clav.volume(),'clavicle_components':len(after_clav.decompose()),
        'full_geometry_verified':False},indent=2))
