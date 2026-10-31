"""M2 end interfaces: wider TPU endplates, Bowden shoulder fittings and C saddles.
Separate candidate to preserve previously audited M1.
"""
from build import *
import os
BASE=OUT/'elbow_connection';DEST=BASE/'muscle_mounts/end_fittings'

def load(path):
    v,f=read_stl(str(path));return mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))

def nearest(s,p):
    v=np.asarray(s.to_mesh64().vert_properties)[:,:3];return v[np.argmin(np.linalg.norm(v-p,axis=1))]

def run():
    DEST.mkdir(exist_ok=True);hum=load(BASE/'humerus_bored.stl');fore=load(BASE/'forearm_bridged.stl');mds=json.loads((OUT/'muscle_parameters.json').read_text(encoding='utf8'));parts={};owners={};recipes=[]
    for md in mds:
        name=md['id'];A=np.array(md['A_sheath_reaction_mm']);C=np.array(md['C_wire_terminal_and_tendon_mm']);B=np.array(md['B_bone_insertion_mm']);z=unit(md['axis']);F=frame(z);x,y=F[:,0],F[:,1];c=A-z*7
        tpu,params=muscle(name,A,C,B,16,12,end_half_width=11,mount_offset=7)
        # A bridge with clearance for M3 heads and a central OD6 ferrule.
        plate=mf.Manifold.cube([26,14,4],True).transform(np.c_[F,c]);links=[]
        for sign in [-1,1]:
            tip=c+x*sign*11;root=nearest(hum,A-z*12+x*sign*6);links.append(rod(root,tip,3,16))
        holes=[rod(c-z*8,c+z*8,2.1,24)]
        for sign in [-1,1]:
            p=c+x*sign*7;holes.append(rod(p-z*8,p+z*8,1.7,24))
        hum=union([hum,plate,*links])-union(holes)
        eye_axis=x;eye=rod(B-x*3,B+x*3,6,24);root=nearest(fore,B);anchor=B+unit(root-B)*4.5;fore=union([fore,eye,rod(root,anchor,3,16)])-rod(B-x*5,B+x*5,1.7,24)
        # Compression shoulder: OD4 tube butts into OD4.2 counterbore;
        # central ID1.4 passage transmits only the steel wire.
        ferrule=union([rod(A-z*17,A-z*9.1,3,24),rod(A-z*9.1,A-z*5.2,2,24)])
        ferrule=ferrule-rod(A-z*19,A-z*3,.7,24)-rod(A-z*18,A-z*13.5,2.1,24)
        saddle=mf.Manifold.cube([22,12,2],True).transform(np.c_[F,C+z*3.6]);eyep=C+y*8+z*6.5
        saddle=union([saddle,rod(eyep-x*1.5,eyep+x*1.5,3.5,24)])-rod(eyep-x*3,eyep+x*3,1.7,24)-rod(C-z*3,C+z*12,.7,24)
        for sign in [-1,1]:
            p=C+x*sign*7;saddle=saddle-rod(p-z*8,p+z*8,1.7,24)
        terminal=rod(C+z*4.7,C+z*8.7,2.2,24)-rod(C+z*3,C+z*10,.7,24)
        p=C+z*6.7;terminal=terminal-rod(p+x*.4,p+x*4,1.05,24)
        pieces={'TPU':tpu,'A_ferrule':ferrule,'C_saddle':saddle,'C_wire_terminal':terminal}
        # Rigid sleeves span the 2.5 mm A-plate/endplate gap without crushing lattice.
        for sign in [-1,1]:
            p=A+x*sign*7;pieces['A_spacer_'+str(sign)]=rod(p-z*4.9,p-z*2.6,2.8,24)-rod(p-z*6,p,1.7,24)
            hum=hum-rod(p-z*5.05,p-z*2.45,3.0,24)
        for suffix,s in pieces.items():parts[name+'_'+suffix]=s;owners[name+'_'+suffix]='humerus' if suffix.startswith('A_') else 'neutral_C' if suffix.startswith('C_') else 'TPU'
        params.update({'endplate_half_width_mm':11,'mount_hole_half_spacing_mm':7,'C_tendon_eye_center_mm':eyep.tolist(),'C_wire_terminal':'metal barrel with longitudinal1.4mm bore and radial nominal M2.5 set-screw pilot; threads/clamp hardware not validated','A_sheath_stop':'PTFE4OD/2ID candidate butts at upstream shoulder; ferrule flange transfers compression to bone plate. Tube pullout restraint still needed when slack.'});recipes.append(params)
    parts.update({'humerus_mounted':hum,'forearm_mounted':fore});owners.update({'humerus_mounted':'humerus','forearm_mounted':'forearm'});cad=[]
    for n,s in parts.items():cad.append({'name':n,**stl(DEST/(n+'.stl'),s)})
    overlaps=[]
    for (an,aa),(bn,bb) in itertools.combinations(parts.items(),2):
        vol=float((aa^bb).volume())
        if vol>.001:overlaps.append({'parts':[an,bn],'volume_mm3':vol})
    dump(DEST/'cad_input.json',{'parts':cad,'units':'mm'});dump(DEST/'interfaces.json',{'recipes':recipes,'owners':owners,'components':{n:len([s for s in p.decompose() if s.volume()>1e-6]) for n,p in parts.items()},'neutral_overlaps':overlaps,'status':'nominal interface candidate; not tested contraction or retention','unmodeled_required_items':['A/C M3 screws/nuts and washer seats','C radial set screw and threads','B/C tendon crimp or clevis','Bowden slack pullout restraint']})
    tree=ET.parse(BASE/'arm.xml');r=tree.getroot()
    for mesh in r.findall('./asset/mesh'):mesh.set('file',os.path.relpath((BASE/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for n,file in [('humerus','humerus_mounted.stl'),('ulna','forearm_mounted.stl'),('elbow_flexor','elbow_flexor_TPU.stl'),('elbow_extensor','elbow_extensor_TPU.stl')]:r.find(f'./asset/mesh[@name="{n}"]').set('file',file)
    origins=json.loads((OUT/'body_frames.json').read_text(encoding='utf8'))['origins_mm']
    for params in recipes:
        p=np.array(params['C_tendon_eye_center_mm'])-origins['humerus'];r.find(f'.//site[@name="{params["id"]}_C"]').set('pos',' '.join(str(v/1000) for v in p))
    for n,s in parts.items():
        if n in ['humerus_mounted','forearm_mounted'] or n.endswith('_TPU'):continue
        ET.SubElement(r.find('asset'),'mesh',name=n,file=n+'.stl',scale='.001 .001 .001');group='4' if owners[n]=='neutral_C' else '3'
        ET.SubElement(r.find('.//body[@name="humerus"]'),'geom',name=n,type='mesh',mesh=n,pos=' '.join(str(-v/1000) for v in origins['humerus']),rgba='.3 .5 .65 1',group=group,contype='0',conaffinity='0',mass='.001')
    ET.indent(r);tree.write(DEST/'arm.xml',encoding='utf8',xml_declaration=True);print('M2 neutral overlaps',overlaps,flush=True)

if __name__=='__main__':run()
