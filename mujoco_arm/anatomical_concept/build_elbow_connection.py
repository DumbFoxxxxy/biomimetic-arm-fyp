"""E1 physical pin/clevis elbow on the independently fitted anatomical arm.
Ulna and radius are joined by explicit bridges as one first-prototype forearm.
Nominal clearances only; no load, fit or real actuation acceptance.
"""
from build_GH_connection import load,clean
from build_SC_connection import unit,rod,write
from repair_inputs import *
from scipy.spatial.transform import Rotation
import os
BASE=OUT/'joint_model/SC_connection/AC_compact/GH_connection'
DEST=BASE/'elbow_connection'

def run():
    DEST.mkdir(exist_ok=True)
    definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'))
    fit=json.loads((OUT/'joint_candidates.json').read_text(encoding='utf8'))
    inventory=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inventory['models']}
    E=np.array(definition['centres_source_mm']['elbow'])-SHIFT;G=np.array(definition['centres_source_mm']['GH'])-SHIFT
    axis=unit(fit['elbow']['axis_source']);up=unit(G-E);up=unit(up-axis*(up@axis));side=np.cross(axis,up)
    hum,hv,hf=load(BASE/'humerus_GH_ball.obj');ulna,uv,uf=load(ROOT/byid['FJ3391']['source'],True);radius,rv,rf=load(ROOT/byid['FJ3349']['source'],True)
    down=unit((uv.mean(0)+rv.mean(0))/2-E)
    def bone_root(vertices,faces,direction,distance):
        # Exact edge/plane crossings avoid a vertex-density-dependent root.
        edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
        a=vertices[edges[:,0]];b=vertices[edges[:,1]];da=(a-E)@direction-distance;db=(b-E)@direction-distance
        mask=(da*db<=0)&(abs(da-db)>1e-10);a=a[mask];b=b[mask];t=da[mask]/(da[mask]-db[mask]);points=a+t[:,None]*(b-a)
        if len(points)<3:raise RuntimeError('No bone cross-section at '+str(distance))
        return np.unique(np.round(points,8),axis=0).mean(0)
    hroot=bone_root(hv,hf,up,43)
    local_vertices=np.concatenate([hv,uv,rv]);near=local_vertices[np.linalg.norm(local_vertices-E,axis=1)<32]
    half=max(21.,float(np.ceil(np.max(abs((near-E)@axis)))+3));ear_width=6.
    def cyl(a,b,r):return rod(E+axis*a,E+axis*b,r)
    def slab(a,b):
        x=unit(up);y=np.cross(axis,x);frame=np.stack([x,y,axis],1)
        return mf.Manifold.cube([1000,1000,b-a]).translate([-500,-500,a]).transform(np.c_[frame,E])
    # Fork side rails stay lateral to the forearm's sweep, joining above the joint.
    fork=hum;rails=[]
    for sign in [-1,1]:
        lo,hi=sorted([sign*half,sign*(half+ear_width)])
        boss=cyl(lo,hi,10)
        start=E+axis*sign*(half+ear_width/2)+up*4
        corner=E+axis*sign*(half+ear_width/2)+up*43
        boss=boss+rod(start,corner,4)+rod(corner,hroot,4)
        boss=boss^slab(-half-ear_width,half+ear_width)
        fork=fork+boss;rails.append([start.tolist(),corner.tolist(),hroot.tolist()])
    # A real central hub and two cross-bridges implement the rigid forearm assumption.
    hub=cyl(-8,8,8);forearm=ulna+radius+hub;bridges=[]
    uroot=bone_root(uv,uf,down,30)
    forearm=forearm+rod(E+down*4,uroot,3.5)
    for distance in [42,155]:
        a=bone_root(uv,uf,down,distance);b=bone_root(rv,rf,down,distance)
        forearm=forearm+rod(a,b,3);bridges.append([a.tolist(),b.tolist()])
    # Local recess clears the rotating hub while retaining source condylar outline.
    fork=fork-cyl(-8.35,8.35,8.35)
    # Clearance of the added hub-to-ulna link along the intended elbow path.
    cuts=[]
    for angle in range(0,151,2):
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix()
        cuts.append(rod(E+R@(down*4),E+R@(uroot-E),3.85))
    fork=fork-mf.Manifold.batch_boolean(cuts,mf.OpType.Add)
    fork=fork-cyl(-half-ear_width-1,half+ear_width+1,5.85)
    bore=cyl(-half-25,half+25,3.3);fork=fork-bore
    forearm=forearm-cyl(-8.1,8.1,4.2)-bore
    bush=cyl(-7.8,7.8,4.)-cyl(-9,9,3.15)
    # Explicit local assembly opening: insert forearm/bush from the side, before
    # sleeves and pin. The real fork rails/ears remain and carry the joint.
    local=forearm^mf.Manifold.sphere(39,32).translate(E)
    offset=mf.Manifold.batch_boolean([local]+[local.translate(d*s*.3) for d in [axis,up,side] for s in [-1,1]],mf.OpType.Add)
    opening=mf.Manifold.batch_boolean([offset.translate(side*distance) for distance in range(46)],mf.OpType.Add)
    opening=opening^mf.Manifold.sphere(39,32).translate(E)
    before_opening=fork.volume();fork=fork-opening;opening_removed=before_opening-fork.volume()
    parts={'humerus_GH_elbow_fork':fork,'forearm_elbow_hub':forearm,'elbow_bearing_bush':bush}
    for name,sign in [('left',-1),('right',1)]:
        a,b=sorted([sign*8.2,sign*(half+ear_width)])
        parts['elbow_spacer_'+name]=cyl(a,b,5.5)-bore
    # Custom shoulder pin envelope: 6mm shank, 10mm head; nominal M6 retaining nut.
    outer=half+ear_width
    pin=cyl(-outer-1.3,outer+9.3,3)+cyl(-outer-5.3,-outer-1.3,5)
    left=cyl(-outer-1.15,-outer-.15,9)-bore
    right=cyl(outer+.15,outer+1.15,9)-bore
    x=unit(up);y=np.cross(axis,x);frame=np.stack([x,y,axis],1)
    nut=mf.Manifold.cylinder(5,10/np.sqrt(3),10/np.sqrt(3),6).transform(np.c_[frame,E+axis*(outer+1.35)])-bore
    parts.update(elbow_shoulder_pin=pin,elbow_washer_left=left,elbow_washer_right=right,elbow_retaining_nut=nut)
    parts={name:clean(s,name) for name,s in parts.items()}
    rows=[write(DEST/(name+'.stl'),s) for name,s in parts.items()];overlaps=[]
    for (an,a),(bn,b) in itertools.combinations(parts.items(),2):
        volume=max(0,float((a^b).volume()))
        if volume>.001:overlaps.append({'parts':[an,bn],'volume_mm3':volume})
    owners={name:'elbow' if name in ['forearm_elbow_hub','elbow_bearing_bush'] else 'GH' for name in parts}
    aliases={'humerus_GH_elbow_fork':'FJ3368','forearm_elbow_hub':'FJ3391'}
    original_forearm=ulna+radius
    report={'scope':'E1 real clevis/shoulder-pin candidate; not fitted, load-tested or driven','centre_source_mm':(E+SHIFT).tolist(),'axis_source':axis.tolist(),'fork_half_gap_mm':half,'ear_width_mm':ear_width,'hub_width_mm':16,'hub_outer_radius_mm':8,'hub_recess_radial_clearance_mm':.35,'bush_outer_radius_mm':4,'bush_inner_radius_mm':3.15,'pin_radius_mm':3,'nominal_axial_hub_to_spacer_gap_mm':.2,'custom_pin_shank_length_mm':2*outer+9.6,'fastener':'6mm shoulder pin envelope, nominal M6 nut; thread not modeled','humerus_rails_display_mm':rails,'forearm_bridges_display_mm':bridges,'components':{r['name']:r['components'] for r in rows},'part_owners':owners,'bone_aliases':aliases,'replaced_source_geoms':['FJ3368','FJ3391','FJ3349'],'removed_geoms':['FJ3349'],'forearm_policy':'Ulna and radius are geometrically joined by two actual 6mm diameter bridges and represented by FJ3391. No pronation/supination in this first assembly. FJ3349 is NOT deleted material: its entire source solid is present inside forearm_elbow_hub except recorded bore cuts.','neutral_overlaps':overlaps,'geometry_delta':{'humerus_removed_mm3':float((hum-parts['humerus_GH_elbow_fork']).volume()),'humerus_added_mm3':float((parts['humerus_GH_elbow_fork']-hum).volume()),'forearm_removed_mm3':float((original_forearm-parts['forearm_elbow_hub']).volume()),'forearm_added_mm3':float((parts['forearm_elbow_hub']-original_forearm).volume())},'surface_simplification_tolerance_mm':.01,'frozen':False,'pending':['whole-arm solid motion','actual axial retention and hardware assembly','muscle attachment and tendon clearance','manufacturing tolerance, strength and loads']}
    report.update(custom_pin_shank_length_mm=2*outer+10.6,fork_sleeve_bore_radius_mm=5.85,sleeve_outer_radius_mm=5.5,washer_outer_radius_mm=9,washer_thickness_mm=1,bone_root_method='triangle-edge / plane crossings at recorded axial stations',assembly_opening={'construction':'local forearm swept translation, physical contour relief','direction':side.tolist(),'distance_mm':45,'construction_step_mm':1,'local_radius_mm':39,'offset_axis_translations_mm':.3,'humeral_material_removed_mm3':float(opening_removed)},hub_link_clearance_sweep={'radius_mm':3.85,'angle_deg':[0,150],'construction_step_deg':2})
    (DEST/'assembly.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    (DEST/'cad_input.json').write_text(json.dumps({'parts':rows,'units':'mm'}),encoding='utf8')
    tree=ET.parse(BASE/'arm.xml');root=tree.getroot();root.find('compiler').set('meshdir','.')
    for mesh in root.findall('./asset/mesh'):mesh.set('file',os.path.relpath((BASE/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for name,fid in aliases.items():root.find(f'./asset/mesh[@name="{fid}"]').set('file',name+'.stl')
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='geom' and child.get('name')=='FJ3349':parent.remove(child)
    for name in parts:
        if name in aliases:continue
        ET.SubElement(root.find('asset'),'mesh',name=name,file=name+'.stl',scale='.001 .001 .001')
        owner=owners[name];centre=E if owner=='elbow' else G
        ET.SubElement(root.find(f'.//body[@name="{owner}"]'),'geom',name=name,type='mesh',mesh=name,pos=' '.join(str(x/1000) for x in -centre),rgba='.3 .55 .65 1',group='3',contype='0',conaffinity='0',mass='.001')
    ET.indent(root);tree.write(DEST/'arm.xml',encoding='utf8',xml_declaration=True)
    print('elbow half gap',half,'parts',report['components'],'overlaps',overlaps,flush=True)

if __name__=='__main__':run()
