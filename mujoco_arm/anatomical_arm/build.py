"""A1 source-size anatomical prototype. Millimetres unless explicitly noted.

Run with D:/MujocoEnv/Scripts/python.exe. No source file is modified.
"""
from pathlib import Path
import sys, json, hashlib, struct, itertools, csv
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'骨骼/_optimize'),str(ROOT/'tmp/shoulder_motion_design_20261029/vendor')]
from decimate import read_stl,read_obj
import manifold3d as mf
OUT=ROOT/'mujoco_arm/exports/anatomical_arm/A1'
NAMES={'clavicle':'锁骨','scapula':'肩胛骨','humerus':'肱骨','ulna':'尺骨','radius':'桡骨'}

def dump(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8')

def unit(v):
    return np.asarray(v)/np.linalg.norm(v)

def stl(path,solid):
    mesh=solid.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3];f=np.asarray(mesh.tri_verts)
    tri=v[f];n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.maximum(np.linalg.norm(n,axis=1),1e-20)[:,None]
    dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);a=np.zeros(len(f),dtype=dt);a['n']=n;a['v']=tri
    path.write_bytes(b'A1 geometry; millimetres'.ljust(80,b' ')+struct.pack('<I',len(f))+a.tobytes())
    return {'vertices':v.tolist(),'faces':f.tolist(),'volume_mm3':solid.volume(),'components':len(solid.decompose())}

def rod(a,b,r,n=12):
    a=np.asarray(a);b=np.asarray(b);z=unit(b-a);x=unit(np.cross(z,[1,0,0] if abs(z[0])<.8 else [0,1,0]));y=np.cross(z,x)
    return mf.Manifold.cylinder(float(np.linalg.norm(b-a)),r,r,n).transform(np.c_[np.stack([x,y,z],axis=1),a])

def ball(p,r):return mf.Manifold.sphere(r,24).translate(p)
def union(items):return mf.Manifold.batch_boolean(items,mf.OpType.Add)
def frame(z):
    z=unit(z);x=unit(np.cross(z,[0,1,0]));return np.stack([x,np.cross(z,x),z],axis=1)

def fit_sphere(p):
    c=np.linalg.lstsq(np.c_[2*p,np.ones(len(p))],(p*p).sum(1),rcond=None)[0];r=np.sqrt(c[3]+c[:3]@c[:3])
    sol=least_squares(lambda a:np.linalg.norm(p-a[:3],axis=1)-a[3],np.r_[c[:3],r],loss='soft_l1')
    err=np.linalg.norm(p-sol.x[:3],axis=1)-sol.x[3]
    return sol.x[:3],float(sol.x[3]),float(np.sqrt(np.mean(err**2)))

def audit_source():
    bones={};rows=[];raw={}
    for en,cn in NAMES.items():
        p=next((ROOT/'骨骼').glob('part - '+cn+'-*.STL'));v,f=read_stl(str(p));ov,of=read_obj(str(ROOT/'骨骼'/(cn+'.obj')))
        shift=v.mean(0)-ov.mean(0) # means can differ from unwelded seams; bbox translation is invariant
        shift=v.min(0)-ov.min(0);bboxerr=np.max(np.abs((ov.max(0)+shift)-v.max(0)))
        # Original OBJ has unwelded seams. Weld using same utility as STL.
        from decimate import weld_stl
        ov,of=weld_stl(ov[of]);s=mf.Manifold(mf.Mesh(ov.astype(np.float32),of.astype(np.uint32)))
        if s.status()!=mf.Error.NoError:raise RuntimeError(f'{en}: source mesh not closed: {s.status()}')
        s=s.translate(shift);raw[en]=np.asarray(s.to_mesh().vert_properties)[:,:3].astype(float)
        bones[en]=s
        rows.append({'bone':en,'source':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'obj_sha256':hashlib.sha256((ROOT/'骨骼'/(cn+'.obj')).read_bytes()).hexdigest(),'native_cad_sha256':hashlib.sha256((ROOT/'骨骼'/(cn+'.SLDPRT')).read_bytes()).hexdigest(),'triangles':len(f),'obj_to_stl_translation_mm':shift.tolist(),'bbox_difference_mm':float(bboxerr),'closed_components':len(s.decompose()),'native_cad':str((ROOT/'骨骼'/(cn+'.SLDPRT')).relative_to(ROOT))})
    # Explicit surface selection in source coordinates, retained for reproducibility.
    h=raw['humerus'];sel=(h[:,2]>1300)&(h[:,0]>105);gh,r,rms=fit_sphere(h[sel])
    cl=raw['clavicle'];sc=cl[cl[:,0]>np.percentile(cl[:,0],96)].mean(0);ac=cl[cl[:,0]<np.percentile(cl[:,0],4)].mean(0)
    distal=h[h[:,2]<1054];e=distal.mean(0)
    # A fixed transverse axis is fitted to the distal articular band by PCA.
    _,_,vt=np.linalg.svd(distal-distal.mean(0),full_matrices=False);axis=vt[0];axis*=np.sign(axis[0])
    T=np.diag([-1.,-1.,1.]);trans=-T@sc
    conv=lambda p:T@p+trans
    landmarks={k:conv(v).tolist() for k,v in {'SC':sc,'AC':ac,'GH':gh,'E':e}.items()}
    landmarks['elbow_axis']=list(T@axis*-1) # positive flexion towards anterior +Y
    bones={k:s.transform(np.c_[T,trans]) for k,s in bones.items()}
    report={'units':'mm; STL size corroborated by existing SolidWorks metre bounding-box logs','scale':1,'laterality':'engineering lateral direction inferred from SC->GH; source medical laterality unverified','source_to_design_4x4':np.vstack([np.c_[T,trans],[0,0,0,1]]).tolist(),'bones':rows,'landmarks_mm':landmarks,'GH_fit':{'selection':'source Z>1300 and X>105','points':int(sel.sum()),'radius_mm':r,'radial_rms_mm':rms,'uncertainty':'selection-dependent; cartilage and clinical landmarks unavailable'},'SC_AC_fit':'centroids of clavicle 4% end bands; engineering reference, not verified clinical centres','elbow_fit':{'selection':'humerus source Z<1054','method':'centroid and first PCA direction of distal band','points':len(distal),'uncertainty':'condyle-band estimate; articular axis requires visual review'},'original_bone_shape_deviation_mm':0,'humerus_GH_to_elbow_mm':float(np.linalg.norm(gh-e))}
    return bones,report

def muscle(name,A,C,B,width,depth,end_half_width=None,mount_offset=None):
    """Warped BCC grid + diagonal outer skin + continuous bored guide + end blocks."""
    A=np.asarray(A);C=np.asarray(C);z=unit(C-A);F=frame(z);L=np.linalg.norm(C-A);pieces=[];curves=[]
    end_half_width=width*.43 if end_half_width is None else end_half_width
    mount_offset=width*.25 if mount_offset is None else mount_offset
    layers=max(4,int(np.ceil(L/15)));bore=1.2;nodes={}
    def point(t,x,y):
        # Square-to-disc mapping avoids rectangular cross-section envelopes.
        taper=.38+.62*np.sin(np.pi*t)**.75;bow=F[:,1]*5*np.sin(np.pi*t)
        disk=np.array([x*np.sqrt(1-y*y/2),y*np.sqrt(1-x*x/2)])
        return A+(C-A)*t+bow+F[:,0]*disk[0]*width*taper+F[:,1]*disk[1]*depth*taper
    for i in range(layers+1):
        for j,k in itertools.product(range(3),repeat=2):nodes[i,j,k]=point(i/layers,j-1,k-1)
    for i in range(layers):
        for j,k in itertools.product(range(2),repeat=2):
            mid=point((i+.5)/layers,j-.5,k-.5)
            for di,dj,dk in itertools.product(range(2),repeat=3):pieces.append(rod(mid,nodes[i+di,j+dj,k+dk],1.0,8))
    # Outer diagonal mesh, joined to end rings and BCC boundary.
    for i in range(layers):
        for j in range(8):
            theta=j*np.pi/4;theta2=(j+1)*np.pi/4
            p=point(i/layers,np.cos(theta),np.sin(theta));q=point((i+1)/layers,np.cos(theta2),np.sin(theta2));pieces.append(rod(p,q,1.0,8))
            p2=point(i/layers,np.cos(theta2),np.sin(theta2));q2=point((i+1)/layers,np.cos(theta),np.sin(theta));pieces.append(rod(p2,q2,1.0,8))
    channel=[point(t,0,0) for t in np.linspace(0,1,layers+1)]
    for a,b in zip(channel[:-1],channel[1:]):pieces.append(rod(a,b,2.5,16))
    # Thick elliptical end plates integral to lattice. Separate clamp holes.
    for t in (0.,1.):
        centre=point(t,0,0);disc=mf.Manifold.cylinder(5,end_half_width,end_half_width,32,True).scale([1,depth/width,1]).transform(np.c_[F,centre]);pieces.append(disc)
    solid=union(pieces)
    holes=[rod(channel[i]-z*3,channel[i+1]+z*3,bore,16) for i in range(layers)]
    for t in (0.,1.):
        for sign in (-1,1):
            p=point(t,0,0)+F[:,0]*sign*mount_offset;holes.append(rod(p-z*6,p+z*6,1.7,16))
    solid=solid-union(holes)
    # Boolean tangencies can enclose microscopic cavities. Fill only negative
    # closed shells <0.1 mm3; keep the continuous through-bore intact.
    shells=solid.decompose();micro=[s for s in shells if -.1<s.volume()<0]
    if micro:solid=union([s for s in shells if s.volume()>0])
    solid=solid.simplify(.01) # printable/CAD-compatible facets; <=0.01 mm surface change
    post=solid.decompose()
    if any(-.1<s.volume()<0 for s in post):solid=union([s for s in post if s.volume()>0])
    return solid,{'id':name,'A_sheath_reaction_mm':A.tolist(),'C_wire_terminal_and_tendon_mm':C.tolist(),'B_bone_insertion_mm':list(B),'channel_centerline_mm':[x.tolist() for x in channel],'axis':z.tolist(),'active_length_mm':float(L),'half_width_mm':width,'half_depth_mm':depth,'cell_layers':layers,'strut_diameter_mm':2,'wire_bore_diameter_mm':2.4,'mount_hole_diameter_mm':3.4,'wire_diameter_candidate_mm':1,'components':len(solid.decompose()),'filled_numerical_cavities_mm3':-sum(s.volume() for s in micro),'status':'manufacturing specimen; no deformation or actuation validation','force_path':'motor -> Bowden fixed at A -> internal steel wire -> C -> short passive tendon -> B; TPU A-C reacts to A clamp','C_constraint':'no imposed slider; actual bend/translation must be measured','printing':'TPU95A candidate; orient long axis upwards; outer overhangs and end plates require slicer review/support'}

def xml_model(bones,lm,muscles,hardware):
    root=ET.Element('mujoco',model='A1_anatomical_candidate');ET.SubElement(root,'compiler',angle='degree',meshdir='assets');ET.SubElement(root,'option',gravity='0 0 0');ET.SubElement(root,'visual');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1')
    SC=np.array(lm['SC']);AC=np.array(lm['AC']);GH=np.array(lm['GH']);E=np.array(lm['E']);bodies={};origins={}
    def body(name,parent,p):
        origin=origins.get(parent,np.zeros(3));x=ET.SubElement(bodies.get(parent,world),'body',name=name,pos=' '.join(map(str,(p-origin)/1000)));bodies[name]=x;origins[name]=p;return x
    def joint(b,n,axis,ran):ET.SubElement(b,'joint',name=n,axis=' '.join(map(str,axis)),range=ran,damping='1')
    b=body('frame',None,SC);ET.SubElement(b,'geom',name='frame_mount',type='box',pos='-.04 0 -.03',size='.015 .035 .055',rgba='.3 .4 .5 1',contype='0',conaffinity='0')
    b=body('clavicle','frame',SC);joint(b,'sc_elev',[0,-1,0],'-30 60');joint(b,'sc_protract',[0,0,1],'-30 30')
    b=body('scapula','clavicle',AC);joint(b,'ac_up',[0,-1,0],'-20 60')
    b=body('humerus','scapula',GH);joint(b,'gh_flex',[1,0,0],'-40 120');joint(b,'gh_abd',[0,-1,0],'-20 120');joint(b,'gh_rot',unit(E-GH),'-90 70')
    b=body('forearm','humerus',E);joint(b,'elbow',lm['elbow_axis'],'0 145')
    owners={'radius':'forearm','ulna':'forearm','clavicle':'clavicle','scapula':'scapula','humerus':'humerus'}
    for name in bones:
        owner=owners[name];ET.SubElement(asset,'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(bodies[owner],'geom',name=name,type='mesh',mesh=name,pos=' '.join(map(str,-origins[owner]/1000)),rgba='.88 .86 .78 1',group='0',contype='0',conaffinity='0',mass='.1')
    for name,solid,owner in hardware:
        ET.SubElement(asset,'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(bodies[owner],'geom',name=name,type='mesh',mesh=name,pos=' '.join(map(str,-origins[owner]/1000)),rgba='.25 .5 .65 .85',group='3',contype='0',conaffinity='0',mass='.01')
    for md in muscles:
        name=md['id'];ET.SubElement(asset,'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(bodies['humerus'],'geom',name=name,type='mesh',mesh=name,pos=' '.join(map(str,-GH/1000)),rgba='.8 .08 .07 1',group='1',contype='0',conaffinity='0',mass='.01')
        for role,owner in [('A_sheath_reaction_mm','humerus'),('C_wire_terminal_and_tendon_mm','humerus'),('B_bone_insertion_mm','forearm')]:
            p=np.asarray(md[role]);ET.SubElement(bodies[owner],'site',name=name+'_'+role[0],pos=' '.join(map(str,(p-origins[owner])/1000)),size='.003',rgba='1 .7 .1 1')
    for item in json.loads((OUT/'whole_arm_layout.json').read_text(encoding='utf8')):
        name=item['id']+'_VOLUME_ONLY';owner=item['origin_owner'];ET.SubElement(asset,'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(bodies[owner],'geom',name=name,type='mesh',mesh=name,pos=' '.join(map(str,-origins[owner]/1000)),rgba='.7 .3 .15 .35',group='2',contype='0',conaffinity='0',mass='.001')
    tendons=ET.SubElement(root,'tendon')
    for md in muscles:
        t=ET.SubElement(tendons,'spatial',name=md['id']+'_tendon_reference',width='.0005',rgba='.12 .12 .12 1');ET.SubElement(t,'site',site=md['id']+'_C');ET.SubElement(t,'site',site=md['id']+'_B')
    ET.indent(root);ET.ElementTree(root).write(OUT/'arm.xml',encoding='utf-8',xml_declaration=True)
    dump(OUT/'body_frames.json',{'origins_mm':{k:v.tolist() for k,v in origins.items()},'bone_owners':owners,'visualization':'bone-only prescribed kinematics; muscle mesh held in neutral shape on humerus; tendon length visual only; no soft simulation','collision':'MuJoCo mesh convex hulls disabled; exact original meshes audited separately'})

def run():
    (OUT/'assets').mkdir(parents=True,exist_ok=True);bones,report=audit_source()
    calibration=OUT/'elbow_calibration.json'
    if calibration.exists():
        cal=json.loads(calibration.read_text(encoding='utf8'));report['initial_landmarks_mm']=dict(report['landmarks_mm']);report['landmarks_mm']['E']=cal['E_mm'];report['landmarks_mm']['elbow_axis']=cal['axis'];report['engineering_axis_calibration']=cal
    dump(OUT/'source_audit.json',report);lm=report['landmarks_mm'];SC,AC,GH,E=[np.array(lm[k]) for k in ['SC','AC','GH','E']];cad=[]
    for n,s in bones.items():
        row=stl(OUT/'assets'/(n+'.stl'),s);cad.append({'name':n,**row});print(n,'components',row['components'],flush=True)
    axis=np.asarray(lm['elbow_axis']);z=unit(E-GH);side=np.array([0.,1.,0.]);muscles=[]
    recipe_path=Path(__file__).parent/'muscle_recipe.json'
    recipe=json.loads(recipe_path.read_text(encoding='utf8')) if recipe_path.exists() else {'half_width_mm':16,'half_depth_mm':12,'proximal_offset_from_GH_mm':45,'distal_offset_before_elbow_mm':42,'front_back_offset_mm':34,'insertion_distal_to_elbow_mm':45,'insertion_front_back_offset_mm':15}
    for name,sign in [('elbow_flexor',1),('elbow_extensor',-1)]:
        A=GH+z*recipe['proximal_offset_from_GH_mm']+side*sign*recipe['front_back_offset_mm'];C=E-z*recipe['distal_offset_before_elbow_mm']+side*sign*recipe['front_back_offset_mm'];B=E+z*recipe['insertion_distal_to_elbow_mm']+side*sign*recipe['insertion_front_back_offset_mm']
        s,md=muscle(name,A,C,B,recipe['half_width_mm'],recipe['half_depth_mm']);row=stl(OUT/'assets'/(name+'.stl'),s);cad.append({'name':name,**row});muscles.append(md);print(name,'components',row['components'],flush=True)
        # Print-oriented copy is physically identical, with axis pointing +Z.
        A0=np.asarray(md['A_sheath_reaction_mm']);F=frame(np.asarray(md['axis']));local=s.transform(np.c_[F.T,-F.T@A0]);bb=np.array(local.bounding_box());local=local.translate([-bb[0],-bb[1],-bb[2]])
        stl(OUT/'assets'/(name+'_print_mm.stl'),local)
    dump(OUT/'muscle_parameters.json',muscles)
    # Editable nominal hardware references, not certified installed mechanisms.
    hardware=[]
    for name,c,a,owner in [('SC_pin',SC,[0,0,1],'frame'),('AC_pin',AC,[0,-1,0],'clavicle'),('elbow_pin',E,axis,'humerus')]:
        a=np.asarray(a);s=rod(c-a*23,c+a*23,3,24);hardware.append((name,s,owner))
    hardware.append(('GH_ball_reference',ball(GH,10),'humerus'))
    for name,c,a,owner in [('elbow_bearing_reference',E,axis,'forearm'),('AC_bearing_reference',AC,np.array([0,-1,0]),'scapula')]:
        s=rod(c-a*15,c+a*15,7,24)-rod(c-a*18,c+a*18,3.3,24);hardware.append((name,s,owner))
    for n,s,o in hardware:
        row=stl(OUT/'assets'/(n+'.stl'),s);cad.append({'name':n,**row})
    dump(OUT/'hardware.json',{'status':'reference geometry only; no bone pockets, clamp bridges or verified retention yet','parts':[{'id':n,'owner':o,'status':'NOT assembly-ready'} for n,s,o in hardware]})
    layout=[]
    for name,owner,origin,insertion,shape in [('anterior_deltoid','scapula',AC,GH+z*100,'curved sheet'),('posterior_deltoid','scapula',AC+[-30,-30,0],GH+z*100,'curved sheet'),('cuff_internal','scapula',GH+[-35,20,-30],GH+z*45,'curved band'),('cuff_external','scapula',GH+[-35,-20,-30],GH+z*45,'curved band'),('scapular_up','frame',SC+[0,0,40],AC+[-20,0,-50],'fan'),('scapular_down','frame',SC+[0,0,-150],AC+[-30,-10,-40],'fan'),('scapular_protract','frame',SC+[0,40,-80],AC+[-50,0,-70],'sheet'),('scapular_retract','frame',SC+[-40,-30,-40],AC+[-50,-15,-70],'sheet')]:
        layout.append({'id':name,'origin_owner':owner,'origin_candidate_mm':list(origin),'insertion_candidate_mm':list(insertion),'shape':shape,'status':'functional volume/path candidate; attachments not surface-projected or validated'})
    for item in layout:
        a=np.asarray(item['origin_candidate_mm']);b=np.asarray(item['insertion_candidate_mm']);F=frame(unit(b-a));s=mf.Manifold.sphere(1,24).scale([17,8,np.linalg.norm(b-a)/2]).transform(np.c_[F,(a+b)/2]);stl(OUT/'assets'/(item['id']+'_VOLUME_ONLY.stl'),s)
    dump(OUT/'whole_arm_layout.json',layout);dump(OUT/'cad_input.json',{'parts':cad,'units':'mm','format':'faceted bone/lattice geometry; reference hardware'})
    xml_model(bones,lm,muscles,hardware)
    dump(OUT/'manifest.json',{'version':'A1','frozen':False,'replaces_V2S1':False,'status':'work-in-progress anatomical and muscle manufacturing candidate','source_audit':report['bones'],'bone_shape_changes':'none','physical_validation':False})
    print('OUTPUT',OUT,flush=True)

if __name__=='__main__':run()
