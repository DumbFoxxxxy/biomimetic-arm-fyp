"""User-authorized concentric GH mechanism with separately cabled axes.

All three axes still intersect at GH. An 80 mm lateral prototype relocation
provides space for a medial-only flexion bearing and an open lateral yoke.
Concentric supports, groove-separated sector
drums and spatial cable routes are explicit geometry, not joint motors.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import copy,json,sys,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.spatial.transform import Rotation
from shoulder_model import ROOT,SOURCE,EXPORT,v8,add_connector,inputs
from audit_arm_rom import build_arm
from search_gh_anchors import A,B,C
from design_axis_spool import BASE,build_stable_spool
sys.path[:0]=[str(SOURCE/'vendor'),str(ROOT.parent/'cable_design_v7'),str(ROOT.parent)]
import manifold3d as md
from elbow_spool_v7 import build_elbow_model
from stl_io import write_binary_stl

GH_SHIFT=80.*B
AXES=[('GF','gh_flex','scapula','gh_flex',B,-C,-55.,-40.,120.),
      ('GA','gh_abd','gh_flex','gh_abd',C,B,-55.,-20.,120.),
      ('GR','gh_rot','gh_abd','humerus',A,B,15.,-90.,70.),
      ('E','elbow','humerus','ulna',B,C,32.,0.,145.)]

def matrix_quat(R):
    q=Rotation.from_matrix(R).as_quat();return v8.vec(np.r_[q[3],q[:3]])

def annulus_asset(root,name,outer,inner,width,center=(0.,0.,0.)):
    version='_open_v2' if name=='gh_outer_yoke' else '_ports_v2'
    version+=f'_at{np.dot(center,A):g}'
    path=EXPORT/f'{name}_{outer:g}_{inner:g}_{width:g}{version}.stl';EXPORT.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        solid=md.Manifold.cylinder(width,outer,circular_segments=72,center=True)-md.Manifold.cylinder(width+2.,inner,circular_segments=72,center=True)
        solid=solid.transform(np.column_stack([np.column_stack([B,C,A]),np.array(center)]))
        if name=='gh_outer_yoke':
            for axis,radius,hole in [(B,6.,3.2),(C,7.,4.5)]:
                h1=axis;helper=A;h2=np.cross(helper,h1);h2/=np.linalg.norm(h2);h3=np.cross(h1,h2)
                Q=np.column_stack([h2,h3,h1])
                for sign in ([-1.] if np.dot(axis,B)>.9 else [-1.,1.]):
                    boss=(md.Manifold.cylinder(8.,radius,circular_segments=40,center=True)-md.Manifold.cylinder(10.,hole,circular_segments=40,center=True)).transform(np.column_stack([Q,sign*40*axis]))
                    solid=solid+boss
                bore=md.Manifold.cylinder(100.,hole,circular_segments=40,center=True).transform(np.column_stack([Q,np.zeros(3)]))
                solid=solid-bore
            # Remove the lateral arc: the upper-arm shaft must cross +B at
            # 90 degree abduction. A full ring or lateral bearing blocks it.
            cutter=md.Manifold.cube((100.,160.,160.),center=True).translate([58.,0.,0.])
            solid=solid-cutter.transform(np.column_stack([np.column_stack([B,C,A]),np.zeros(3)]))
        elif name=='gh_middle_yoke':
            # Middle rotor carries its shafts; outer parent has clearance
            # bores at C +/-40 rather than solid pins through this ring.
            Q=np.column_stack([B,A,-C])
            for sign in [-1.,1.]:
                boss=md.Manifold.cylinder(8.,6.,circular_segments=40,center=True).transform(np.column_stack([Q,sign*30*C]))
                solid=solid+boss
        mesh=solid.to_mesh();points=np.asarray(mesh.vert_properties)[:,:3]
        write_binary_stl(path,(points@v8.DESIGN_TO_MJ.T*.001)[np.asarray(mesh.tri_verts)],'Concentric GH ring: body-local MuJoCo metres')
    ET.SubElement(root.find('asset'),'mesh',{'name':name,'file':str(path)})

def transformed_node(node,Q,t,phase,axial):
    item=copy.deepcopy(node);oldpos=np.fromstring(item.get('pos','0 0 0'),sep=' ')
    item.set('pos',v8.vec(Q@(phase@(oldpos+np.array([axial,0.,0.])))+t))
    if item.tag=='geom':
        oldquat=np.fromstring(item.get('quat','1 0 0 0'),sep=' ')
        oldR=Rotation.from_quat(np.r_[oldquat[1:],oldquat[0]]).as_matrix()
        item.set('quat',matrix_quat(Q@phase@oldR))
    return item

def add_axis(root,prefix,joint_name,parent_name,child_name,axis,outward,offset,minimum,maximum):
    parent=root.find(f'.//body[@name="{parent_name}"]');child=root.find(f'.//body[@name="{child_name}"]')
    cfg=v8.load_config();gh=np.array([159.06,187.54,7.30])+GH_SHIFT;ac=np.array(cfg['geometry_mm']['AC']);elbow=np.array([184.34,-101.35,7.30])
    # All parent/child frames coincide at GH, except the scapula frame at AC
    # and the ulna frame at the elbow.
    parent_origin=(gh-ac if parent_name=='scapula' else (elbow-np.array([159.06,187.54,7.30]) if joint_name=='elbow' else np.zeros(3)))
    center_parent=parent_origin+axis*offset;center_child=axis*offset
    # Align the old X hinge axis with the new axis. Negative old Y points
    # outboard, keeping the two inlet guides away from the torso and shaft.
    D=np.column_stack([axis,-outward,-np.cross(axis,outward)])
    Q=v8.DESIGN_TO_MJ@D
    phase=v8.rotation_matrix(np.array([1.,0.,0.]),np.radians(-minimum))
    radius=25. if prefix=='GR' else 18.
    params={**BASE,'radius_mm':radius,'elbow_range_deg':[0.,maximum-minimum]}
    if prefix=='GR':params.update(entry_radius_mm=3.,entry_width_mm=2.,entry_guide_mm=[0.,19.,21.],external_anchor_mm=[0.,52.,13.],sector_side_angle_E01_deg=120.,sector_side_angle_E02_deg=195.)
    source=build_stable_spool(params).getroot()
    if prefix=='GR':
        for rid in ['E01','E02']:
            entry=source.find(f'.//geom[@name="entry_{rid}"]');bridge=source.find(f'.//site[@name="bridge_{rid}"]')
            p=np.fromstring(entry.get('pos'),sep=' ');q=np.fromstring(bridge.get('pos'),sep=' ')
            bridge.set('pos',v8.vec(p+(q-p)*8./25.))
    name_map={}
    for node in source.iter():
        if node.get('name'):name_map[node.get('name')]=f'{prefix}_{node.get("name")}'
    source_parent=source.find('.//body[@name="humerus"]');source_child=source.find('.//body[@name="ulna"]')
    for owner,source_body,translation,rot in [(parent,source_parent,v8.design_to_mj(center_parent),np.eye(3)),(child,source_child,v8.design_to_mj(center_child),phase)]:
        for node in source_body:
            if node.tag not in ['geom','site'] or node.get('name') in ['humerus_link','ulna_link']:continue
            rid='E01' if 'E01' in node.get('name','') else 'E02'
            axial=.003 if rid=='E01' else -.003
            item=transformed_node(node,Q,translation,rot,axial)
            item.set('name',name_map[node.get('name')])
            if item.tag=='geom' and 'sector_' in node.get('name',''):
                item.set('size',f'{radius*.001} .0015')
            if item.tag=='geom':
                item.set('contype','0');item.set('conaffinity','0')
                item.attrib.pop('material',None);item.set('rgba','.6 .7 .85 1' if 'sector_' in node.get('name','') else '.85 .35 .2 1')
                if 'sector_' in node.get('name',''):item.set('rgba','0 0 0 0');item.set('group','3')
            owner.append(item)
    for node in source.find('worldbody'):
        if node.tag!='site' or not node.get('name','').startswith('external_'):continue
        rid='E01' if 'E01' in node.get('name') else 'E02';axial=.003 if rid=='E01' else -.003
        item=transformed_node(node,Q,v8.design_to_mj(center_parent),np.eye(3),axial);item.set('name',name_map[node.get('name')]);parent.append(item)
    # Actual hollow two-groove drum. Invisible cylinders above are only the
    # analytic tendon-wrapping proxies; they must not fill the physical bore.
    path=EXPORT/f'{prefix}_drum_r{radius:g}_d{offset:g}_grooves_v3.stl'
    if not path.exists():
        solid=md.Manifold.cylinder(15.,radius-.5,circular_segments=96,center=True)
        for level in [-5.,0.,5.]:solid=solid+md.Manifold.cylinder(2.,radius+1.,circular_segments=96,center=True).translate([0.,0.,level])
        solid=solid-md.Manifold.cylinder(17.,3.2,circular_segments=40,center=True)
        mesh=solid.to_mesh();uu=outward;vv=np.cross(axis,uu);p=np.asarray(mesh.vert_properties)[:,:3]@np.column_stack([uu,vv,axis]).T+center_child
        write_binary_stl(path,(p@v8.DESIGN_TO_MJ.T*.001)[np.asarray(mesh.tri_verts)],'Hollow two-groove cable drum; wrapping proxies excluded')
    ET.SubElement(root.find('asset'),'mesh',{'name':f'{prefix}_drum','file':str(path)})
    ET.SubElement(child,'geom',{'name':f'{prefix}_drum','type':'mesh','mesh':f'{prefix}_drum','density':'0','contype':'0','conaffinity':'0','rgba':'.45 .5 .6 1'})
    for node in source.find('tendon'):
        item=copy.deepcopy(node);item.set('name',f'{prefix}_{node.get("name")}')
        # The motor-to-parent-guide section is fixed within the carrier and
        # contributes a constant supply length. Start the active tendon at the
        # guide exit; a degenerate doubled-back fixed guide is not a moving
        # transmission and must not create spurious wrap-branch jumps.
        for fixed_supply_node in list(item)[:2]:item.remove(fixed_supply_node)
        for ref in item:
            for attr in ['site','geom','sidesite']:
                if attr in ref.attrib:ref.set(attr,name_map[ref.get(attr)])
        root.find('tendon').append(item)
        rid=item.get('name');ET.SubElement(root.find('actuator'),'position',{'name':f'pull_{rid}','tendon':rid,'kp':'10000','kv':'150','ctrllimited':'true','ctrlrange':'0 2','forcelimited':'true','forcerange':'-200 0'})
    # Dedicated inlet guide mounts; the steel supply side may be led through
    # sheaths whose distal ferrules are fixed to this parent carrier.
    for rid in ['E01','E02']:
        entry=parent.find(f'geom[@name="{prefix}_entry_{rid}"]');center=np.fromstring(entry.get('pos'),sep=' ')
        local_center=v8.mj_to_design(center)
        if prefix=='GR':
            radial=local_center-center_parent;radial-=axis*np.dot(radial,axis);radial/=np.linalg.norm(radial)
            start=15*A+30*radial
            add_connector(parent,f'{prefix}_{rid}_guide_root',-10*A+25*radial,-10*A+30*radial,2.)
            add_connector(parent,f'{prefix}_{rid}_guide_stem',-10*A+30*radial,start,2.)
        else:start=center_parent+25*outward
        add_connector(parent,f'{prefix}_{rid}_guide_mount',start,local_center,2.)
    if prefix in ['GF','GA']:
        direction=1. if offset>0 else -1.
        pivot=parent_origin+axis*(36.*direction);far=center_parent+axis*(15.*direction)
        if prefix=='GF':add_connector(parent,f'{prefix}_fixed_axle',pivot,far,3.)
        else:
            # Attach the fixed inlet carrier around the rotating C shaft.
            add_connector(parent,'GA_bearing_carrier_root',parent_origin-40*C+7*B,parent_origin-40*C+25*B,2.)
            add_connector(parent,'GA_bearing_carrier_side',parent_origin-40*C+25*B,far+25*B,2.)
        add_connector(parent,f'{prefix}_guide_carrier',far,far+25*outward,2.)
        add_connector(parent,f'{prefix}_guide_carrier_return',far+25*outward,center_parent+25*outward,2.)
    elif prefix=='E':
        start=parent_origin-35*A+20*B;via=start+25*C;end=center_parent+25*C
        add_connector(parent,'E_holder_root',start,via,2.)
        add_connector(parent,'E_holder_arm',via,end,2.)
        add_connector(child,'E_rotor_connection',np.zeros(3),center_child,4.)
    return {'prefix':prefix,'joint':joint_name,'parent':parent_name,'child':child_name,'axis_design':axis.tolist(),
            'drum_center_parent_mm':center_parent.tolist(),'drum_center_child_mm':center_child.tolist(),'range_deg':[minimum,maximum],
            'drum_radius_mm':radius,'physical_groove_root_radius_mm':radius-.5,'physical_flange_radius_mm':radius+1.,
            'source_guide_radius_mm':3. if prefix=='GR' else 6.,'supply_interface':'GR: bored fixed Bowden ferrules; other axes: parent bridge coordinates, final sheath clamps still need detailing',
            'routing_scope':'Active bridge-to-drum-to-anchor; upstream motor/sheath segment is not included','groove_plane_offsets_mm':[-3.,3.],
            'parent_origin_mm':parent_origin.tolist(),'basis_old_to_design':D.tolist()}

def build_gimbal():
    tree=build_arm();root=tree.getroot()
    routes=EXPORT/'full_arm_shoulder_routes_candidate.json'
    if routes.exists():
        candidate=json.loads(routes.read_text())
        if candidate['pass']:
            # Redistribute the stabilizing lines when the complete-arm route
            # search selects different eyes. The shoulder-only checkpoint is
            # retained separately; unused old actuators must not remain active.
            for owner in [root.find('tendon'),root.find('actuator')]:
                for node in list(owner):
                    name=node.get('name','').removeprefix('pull_')
                    if name.startswith('S') and name[1:].isdigit():owner.remove(node)
            for owner,stem in [(root.find('worldbody'),'fixed_S'),(root.find('.//body[@name="scapula"]'),'eye_S')]:
                for node in list(owner):
                    if node.tag=='site' and node.get('name','').startswith(stem):owner.remove(node)
            for row in candidate['routes']:
                rid=row['id'];fixed=root.find(f'.//site[@name="fixed_{rid}"]')
                if fixed is None:
                    fixed=ET.SubElement(root.find('worldbody'),'site',{'name':f'fixed_{rid}','size':'.002'})
                    ET.SubElement(root.find('.//body[@name="scapula"]'),'site',{'name':f'eye_{rid}','pos':v8.vec(v8.design_to_mj(row['eye_local_ac_mm'])),'size':'.0015'})
                    cable=ET.SubElement(root.find('tendon'),'spatial',{'name':rid,'width':'.001','rgba':'.95 .4 .1 1'})
                    ET.SubElement(cable,'site',{'site':f'fixed_{rid}'});ET.SubElement(cable,'site',{'site':f'eye_{rid}'})
                    ET.SubElement(root.find('actuator'),'position',{'name':f'pull_{rid}','tendon':rid,'kp':'10000','kv':'150','ctrllimited':'true','ctrlrange':'0 2','forcelimited':'true','forcerange':'-200 0'})
                fixed.set('pos',v8.vec(v8.design_to_mj(row['fixed_eye_world_mm'])))
    gh=root.find('.//body[@name="gh_flex"]');oldpos=np.fromstring(gh.get('pos'),sep=' ')
    gh.set('pos',v8.vec(oldpos+v8.design_to_mj(GH_SHIFT)))
    scap=root.find('.//body[@name="scapula"]');obsolete=scap.find('geom[@name="gh_scapula_web"]')
    if obsolete is not None:scap.remove(obsolete)
    humerus=root.find('.//body[@name="humerus"]');oldball=humerus.find('geom[@name="gh_ball"]')
    if oldball is not None:humerus.remove(oldball)
    annulus_asset(root,'gh_outer_yoke',42.,38.,4.)
    annulus_asset(root,'gh_middle_yoke',32.,28.,4.)
    annulus_asset(root,'gh_axial_bearing',26.,21.,6.,-10*A)
    annulus_asset(root,'gh_axial_adapter',20.5,18.,26.)
    for owner,name,color in [('gh_flex','gh_outer_yoke','.85 .5 .2 1'),('gh_abd','gh_middle_yoke','.2 .75 .75 1'),('gh_abd','gh_axial_bearing','.6 .65 .7 1'),('humerus','gh_axial_adapter','.65 .6 .5 1')]:
        ET.SubElement(root.find(f'.//body[@name="{owner}"]'),'geom',{'name':name,'type':'mesh','mesh':name,'rgba':color,'density':'0','contype':'0','conaffinity':'0'})
    cfg=v8.load_config();oldgh=np.array([159.06,187.54,7.30])-np.array(cfg['geometry_mm']['AC']);center=oldgh+GH_SHIFT
    add_connector(scap,'gh_fixed_medial_mount',oldgh+12*v8.source_basis(cfg)[2],center-70*B,3.)
    for sign in [-1.,1.]:
        add_connector(root.find('.//body[@name="gh_abd"]'),f'gh_axial_support_{sign}',sign*30*B,-10*A+sign*24*B,2.5)
    # Hollow shafts connect the rotating yokes to the outboard drums without
    # occupying the fixed bearing pins. Their bores are actual mesh geometry.
    for name,owner,axis,start,end in [('gf_shaft','gh_flex',B,-55.,-40.),('ga_shaft','gh_abd',C,-55.,-30.)]:
        path=EXPORT/f'{name}_{start:g}_{end:g}_tube.stl'
        if not path.exists():
            length=end-start;helper=A;u=np.cross(helper,axis);u/=np.linalg.norm(u);vv=np.cross(axis,u)
            tube=md.Manifold.cylinder(length,4.,circular_segments=40,center=True)-md.Manifold.cylinder(length+2.,3.2,circular_segments=40,center=True)
            mesh=tube.to_mesh();p=np.asarray(mesh.vert_properties)[:,:3]@np.column_stack([u,vv,axis]).T+axis*(start+end)/2
            write_binary_stl(path,(p@v8.DESIGN_TO_MJ.T*.001)[np.asarray(mesh.tri_verts)],'GH hollow rotor shaft: body-local metres')
        ET.SubElement(root.find('asset'),'mesh',{'name':name,'file':str(path)})
        ET.SubElement(root.find(f'.//body[@name="{owner}"]'),'geom',{'name':name,'type':'mesh','mesh':name,'density':'0','contype':'0','conaffinity':'0','rgba':'.6 .65 .7 1'})
    add_connector(root.find('.//body[@name="gh_abd"]'),'ga_positive_journal',30*C,42*C,4.)
    specs=[add_axis(root,*axis) for axis in AXES]
    from refine_gimbal_geometry import refine
    refine(root)
    from prototype_mass_properties import assign
    assign(root)
    ET.indent(tree,space='  ')
    return tree,specs

if __name__=='__main__':
    tree,specs=build_gimbal();model=_load_project_xml(ET.tostring(tree.getroot(),encoding='unicode'))
    tree.write(EXPORT/'gimbal_arm.xml',encoding='utf-8',xml_declaration=True)
    (EXPORT/'gimbal_axis_layout.json').write_text(json.dumps({'gh_shift_design_mm':GH_SHIFT.tolist(),'gh_center_design_mm':(np.array([159.06,187.54,7.30])+GH_SHIFT).tolist(),'axes':specs,'full_geometry_verified':False},indent=2))
    print({'nq':model.nq,'nv':model.nv,'nu':model.nu,'ntendon':model.ntendon})
