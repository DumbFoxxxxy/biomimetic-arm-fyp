"""G1 local engineered GH bearing on the independent anatomical arm.
Slots are actual neck-clearance openings; retained sectors have their own screws.
No source joint centres move, no source files change, no drive/load claim.
"""
from build_AC_connection import load
from build_SC_connection import unit,rod,write
from repair_inputs import *
from scipy.spatial.transform import Rotation
import os

BASE=OUT/'joint_model/SC_connection/AC_compact'
CLOSED='--closed' in sys.argv
PRESERVE_HEAD='--preserve-head' in sys.argv
DEST=BASE/('GH_closed' if CLOSED else 'GH_preserve_head' if PRESERVE_HEAD else 'GH_connection')

def clean(s,name):
    s=s.as_original().simplify(.01);shells=s.decompose();positive=[a for a in shells if abs(a.volume())>1e-10]
    if len(positive)==1:s=positive[0]
    if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid '+name)
    if len(positive)!=1:raise RuntimeError('Disconnected material '+name+' '+str([{'volume':a.volume(),'bounds':a.bounding_box()} for a in positive]))
    return s

def run():
    DEST.mkdir(exist_ok=True);definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'))
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));byid={r['id']:r for r in inv['models']}
    G=np.array(definition['centres_source_mm']['GH'])-SHIFT;E=np.array(definition['centres_source_mm']['elbow'])-SHIFT
    oldscap,sv,sf=load(BASE/'scapula_AC_socket.obj');oldhum,hv,hf=load(ROOT/byid['FJ3368']['source'],True)
    long=unit(E-G);along=(hv-G)@long;root=hv[(along>38)&(along<46)].mean(0)
    if not np.isfinite(root).all():raise RuntimeError('No humeral shaft root')
    z=unit(root-G);x=unit(np.cross(z,[0,-1,0]));y=np.cross(z,x);F=np.stack([x,y,z],1)
    radius=20. if PRESERVE_HEAD else 16.;inner=radius+.4;outer=radius+3.5;neck=4.5;depth=3.;relief=26. if PRESERVE_HEAD else 32.
    def local(s):return s.transform(np.c_[F,G])
    def box(lo,hi):return local(mf.Manifold.cube(np.array(hi)-lo).translate(lo))
    hum=(oldhum-mf.Manifold.sphere(relief,40).translate(G))+mf.Manifold.sphere(radius,40).translate(G)+rod(G,root,neck)
    humeral_braces=[]
    if PRESERVE_HEAD:
        fragments=[s for s in hum.decompose() if abs(s.volume())>1e-10];primary=max(fragments,key=lambda s:s.volume())
        for fragment in fragments:
            if fragment is primary:continue
            v=np.asarray(fragment.to_mesh64().vert_properties)[:,:3];tip=v[np.argmax(np.linalg.norm(v-G,axis=1))]
            a=unit(tip-G);b=unit(root-G);axis=unit(np.cross(a,b));angle=np.arccos(np.clip(a@b,-1,1))
            arc=[G+Rotation.from_rotvec(axis*t).apply(a)*(outer+12) for t in np.linspace(0,angle,max(2,int(np.ceil(np.degrees(angle)/3))+1))]
            points=[tip-a,*arc,root];brace=mf.Manifold.batch_boolean([rod(p,q,2.5) for p,q in zip(points,points[1:])],mf.OpType.Add)
            hum=hum+brace;humeral_braces.append([p.tolist() for p in points])
    shell=mf.Manifold.sphere(outer,40).translate(G)-mf.Manifold.sphere(inner,40).translate(G)
    flex=unit(np.cross(long,[0,-1,0]));abd=unit(np.cross(long,[-1,0,0]));cuts=[]
    if not CLOSED:
        for axis in [flex,abd]:
            for angle in range(0,161,2):
                direction=Rotation.from_rotvec(axis*np.radians(angle)).apply(z)
                cuts.append(rod(G,G+direction*(outer+10),neck+.8))
        swept=mf.Manifold.batch_boolean(cuts,mf.OpType.Add)
        shell=shell-swept
    rear=shell^box([-60,-60,-60],[60,60,-.1]);cap=shell^box([-60,-60,.1],[60,60,depth])
    if CLOSED:
        sectors=[cap^box([-60,-60,-60],[-.1,60,60]),cap^box([.1,-60,-60],[60,60,60])]
    else:sectors=cap.decompose()
    # Preserve every positive cap component; no small useful piece is discarded.
    sectors=[s for s in sectors if abs(s.volume())>1e-10]
    if not sectors:raise RuntimeError('No retained cap material after neck slots')
    base=oldscap-mf.Manifold.sphere(outer+6.5,40).translate(G)
    if PRESERVE_HEAD:
        AC=np.array(definition['centres_source_mm']['AC'])-SHIFT
        base=base+(oldscap^mf.Manifold.sphere(12,32).translate(AC))
    candidates=sv[(np.linalg.norm(sv-G,axis=1)>28)&(np.linalg.norm(sv-G,axis=1)<60)]
    roots=[]
    for sign in [-1,1]:
        target=G-z*35+x*sign*9;attachment=candidates[np.argmin(np.linalg.norm(candidates-target,axis=1))]
        roots.append(attachment.tolist());rear=rear+rod(G-z*(outer-1.5)+x*sign*3,attachment,3.5)
    base=base+rear;parts={};hardware={};bolt_sites=[]
    for i,sector in enumerate(sectors):
        mesh=sector.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3];radial=(v.mean(0)-G);radial-=z*(radial@z);u=unit(radial)
        c=G+u*(outer+2);bolt_sites.append(c.tolist())
        sector=sector+rod(c+z*.1,c+z*depth,3.7)+rod(G+u*(outer-1)+z*.5,c+z*.5,2.)
        sector=sector^box([-60,-60,.1],[60,60,depth])
        sector=sector-mf.Manifold.sphere(inner,40).translate(G)
        base=base+rod(c-z*4,c-z*.1,3.7)+rod(G+u*(outer-2)-z*2,c-z*2,3)
        bore=rod(c-z*14,c+z*10,1.7);sector=sector-bore;base=base-bore
        head=depth+.75
        hardware[f'GH_bolt_{i}']=rod(c+z*(head-12),c+z*head,1.5)+rod(c+z*head,c+z*(head+3),2.75)
        hardware[f'GH_nut_{i}']=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[F,c-z*7.15])-rod(c-z*9,c-z*3,1.6)
        hardware[f'GH_top_washer_{i}']=rod(c+z*(depth+.15),c+z*(depth+.65),3.5)-rod(c+z*(depth-.5),c+z*(depth+1.5),1.7)
        hardware[f'GH_bottom_washer_{i}']=rod(c-z*4.65,c-z*4.15,3.5)-rod(c-z*5,c-z*3,1.7)
        base=base-(rod(c-z*7.4,c-z*4.5,3.45)+rod(c-z*4.9,c-z*3.95,3.75))
        parts[f'GH_cap_{i}']=sector
    for sector in parts.values():
        base=base-mf.Manifold.batch_boolean([sector]+[sector.translate(d*s*.15) for d in [x,y,z] for s in [-1,1]],mf.OpType.Add)
    # Local zero-pose assembly clearance, recorded; motion still needs full audit.
    hum_clear=mf.Manifold.batch_boolean([hum]+[hum.translate(d*s*.15) for d in [x,y,z] for s in [-1,1]],mf.OpType.Add)
    base=base-hum_clear
    braces=[]
    fragments=[s for s in base.decompose() if abs(s.volume())>1e-10]
    if len(fragments)>1:
        primary=max(fragments,key=lambda s:s.volume())
        for fragment in fragments:
            if fragment is primary:continue
            v=np.asarray(fragment.to_mesh64().vert_properties)[:,:3]
            tip=v[np.argmax(np.linalg.norm(v-G,axis=1))];waypoint=G+unit(tip-G)*(outer+8)
            target=candidates[np.argmin(np.linalg.norm(candidates-waypoint,axis=1))]
            brace=rod(tip,waypoint,3.5)+rod(waypoint,target,3.5)
            brace=brace-hum_clear
            for sector in parts.values():brace=brace-sector
            for accessory in hardware.values():brace=brace-accessory
            base=base+brace;braces.append([tip.tolist(),waypoint.tolist(),target.tolist()])
    context=json.loads((BASE/'assembly.json').read_text(encoding='utf8'));context_reliefs=[]
    for name in context['components']:
        if name=='scapula_AC_socket':continue
        existing=load(BASE/(name+'.obj'))[0];overlap=max(0,float((base^existing).volume()))
        if overlap>.001:
            clearance=mf.Manifold.batch_boolean([existing]+[existing.translate(d*s*.15) for d in [x,y,z] for s in [-1,1]],mf.OpType.Add)
            before=base.volume();base=base-clearance;context_reliefs.append({'part':name,'overlap_before_mm3':overlap,'removed_clearance_mm3':float(before-base.volume())})
    trim_regions=[]
    for fragment in base.decompose():
        bounds=np.array(fragment.bounding_box());span=np.linalg.norm(bounds[3:]-bounds[:3])
        if 1e-10<abs(fragment.volume())<.01 and span<2.:
            centre=(bounds[:3]+bounds[3:])/2;radius_trim=span/2+.05
            # Explicit small contour trim of a newly cut floating sliver.
            # Do not discard a detached useful connector or a source bone island.
            base=base-mf.Manifold.sphere(radius_trim,20).translate(centre)
            trim_regions.append({'centre_display_mm':centre.tolist(),'radius_mm':radius_trim,'newly_cut_sliver_volume_mm3':float(fragment.volume())})
    parts={'humerus_GH_ball':hum,'scapula_AC_GH_sockets':base,**parts,**hardware}
    parts={name:clean(s,name) for name,s in parts.items()}
    rows=[write(DEST/(name+'.stl'),s) for name,s in parts.items()];overlaps=[]
    for (an,a),(bn,b) in itertools.combinations(parts.items(),2):
        vol=max(0,float((a^b).volume()))
        if vol>.001:overlaps.append({'parts':[an,bn],'volume_mm3':vol})
    owners={name:'GH' if name=='humerus_GH_ball' else 'AC' for name in parts}
    report={'scope':'G1 engineered GH geometry, independent anatomical concept, not frozen or strength/drive acceptance','centre_source_mm':(G+SHIFT).tolist(),'ball_radius_mm':radius,'socket_inner_radius_mm':inner,'socket_outer_radius_mm':outer,'neck_radius_mm':neck,'neck_axis':z.tolist(),'humeral_root_display_mm':root.tolist(),'scapular_roots_display_mm':roots,'cap_depth_mm':depth,'cap_count':len(sectors),'neck_slot':None if CLOSED else {'construction':'solid swept-cylinder difference','radius_mm':neck+.8,'sweep_step_deg':2,'angles_deg':[0,160],'axes':[flex.tolist(),abd.tolist()],'interpretation':'clearance construction only; not path validation or a general escape proof'},'bolt_sites_display_mm':bolt_sites,'fasteners':'one nominal M3x12, nut and two washers per cap; unthreaded envelopes','components':{r['name']:r['components'] for r in rows},'part_owners':owners,'bone_aliases':{'humerus_GH_ball':'FJ3368','scapula_AC_GH_sockets':'FJ3384'},'neutral_overlaps':overlaps,'geometry_delta':{'humerus_removed_mm3':float((oldhum-parts['humerus_GH_ball']).volume()),'humerus_added_mm3':float((parts['humerus_GH_ball']-oldhum).volume()),'scapula_removed_mm3':float((oldscap-parts['scapula_AC_GH_sockets']).volume()),'scapula_added_mm3':float((parts['scapula_AC_GH_sockets']-oldscap).volume())},'surface_simplification_tolerance_mm':.01,'frozen':False,'pending':['capture directions and escape paths','full solid motion sweep','muscle attachment remapping','assembly access','printed clearance, strength, loads and drive']}
    report['exterior_braces_display_mm']=braces
    report['humeral_tuberosity_braces_display_mm']=humeral_braces
    report['existing_AC_context_reliefs']=context_reliefs
    report['local_clearance_sliver_trims']=trim_regions
    (DEST/'assembly.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    (DEST/'cad_input.json').write_text(json.dumps({'parts':rows,'units':'mm'}),encoding='utf8')
    tree=ET.parse(BASE/'arm.xml');r=tree.getroot();r.find('compiler').set('meshdir','.')
    for mesh in r.findall('./asset/mesh'):mesh.set('file',os.path.relpath((BASE/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for name,fid in report['bone_aliases'].items():r.find(f'./asset/mesh[@name="{fid}"]').set('file',name+'.stl')
    A=np.array(definition['centres_source_mm']['AC'])-SHIFT
    for name in parts:
        if name in report['bone_aliases']:continue
        ET.SubElement(r.find('asset'),'mesh',name=name,file=name+'.stl',scale='.001 .001 .001')
        ET.SubElement(r.find('.//body[@name="AC"]'),'geom',name=name,type='mesh',mesh=name,pos=' '.join(str(v/1000) for v in -A),rgba='.7 .42 .2 1',group='3',contype='0',conaffinity='0',mass='.001')
    ET.indent(r);tree.write(DEST/'arm.xml',encoding='utf-8',xml_declaration=True)
    print('GH caps',len(sectors),'parts',report['components'],'overlaps',overlaps,flush=True)

if __name__=='__main__':run()
