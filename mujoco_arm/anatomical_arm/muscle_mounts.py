"""Create explicit A reaction plates and B tendon eyes on E1 bone parts.
No C slider or fictitious belly motion is introduced.
"""
from build import *
import os
BASE=OUT/'elbow_connection';DEST=BASE/'muscle_mounts'

def load(path):
    v,f=read_stl(str(path));return mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))

def near(s,p):
    v=np.asarray(s.to_mesh64().vert_properties)[:,:3];return v[np.argmin(np.linalg.norm(v-p,axis=1))]

def run():
    DEST.mkdir(exist_ok=True);hum=load(BASE/'humerus_bored.stl');fore=load(BASE/'forearm_bridged.stl');mds=json.loads((OUT/'muscle_parameters.json').read_text(encoding='utf8'));interfaces=[]
    for md in mds:
        A=np.array(md['A_sheath_reaction_mm']);B=np.array(md['B_bone_insertion_mm']);z=unit(md['axis']);F=frame(z);c=A-z*7
        plate=mf.Manifold.cube([18,12,4],True).transform(np.c_[F,c]);supports=[]
        for sign in [-1,1]:
            tip=c+F[:,0]*sign*7;root=near(hum,A-z*12+F[:,0]*sign*6);supports.append(rod(root,tip,3,16))
        # PTFE entrance and two M3 fastener passages aligned with TPU end holes.
        holes=[rod(c-z*8,c+z*8,2.1,24)]
        for sign in [-1,1]:
            p=c+F[:,0]*sign*4;holes.append(rod(p-z*8,p+z*8,1.7,24))
        hum=union([hum,plate,*supports])-union(holes)
        eye_axis=F[:,0];eye=rod(B-eye_axis*3,B+eye_axis*3,6,24);root=near(fore,B);direction=unit(root-B);anchor=B+direction*4.5
        fore=union([fore,eye,rod(root,anchor,3,16)])-rod(B-eye_axis*5,B+eye_axis*5,1.7,24)
        interfaces.append({'muscle':md['id'],'A_reaction_plate_center_mm':c.tolist(),'A_mount_hole_centers_mm':[(A+F[:,0]*sign*4).tolist() for sign in [-1,1]],'A_sheath_entrance_diameter_mm':4.2,'B_eye_center_mm':B.tolist(),'B_eye_axis':eye_axis.tolist(),'B_eye_hole_diameter_mm':3.4,'C':'free TPU distal end; connect active wire terminal and separate passive tendon; no artificial rigid slide'})
    solids={'humerus_mounted':hum,'forearm_mounted':fore};cad=[]
    for n,s in solids.items():cad.append({'name':n,**stl(DEST/(n+'.stl'),s)})
    neutral=[]
    for md in mds:
        s=load(OUT/'assets'/(md['id']+'.stl'))
        for n,bone in solids.items():
            vol=float((s^bone).volume());neutral.append({'muscle':md['id'],'part':n,'overlap_mm3':vol,'pass':vol<=.001})
    dump(DEST/'cad_input.json',{'parts':cad,'units':'mm'});dump(DEST/'interfaces.json',{'interfaces':interfaces,'components':{n:len([x for x in s.decompose() if x.volume()>1e-6]) for n,s in solids.items()},'neutral_TPU_checks':neutral,'qualification':'nominal mount geometry only; Bowden flange retention/clamp and wire end fittings still require mechanical review; no force or deformation acceptance'})
    tree=ET.parse(BASE/'arm.xml');r=tree.getroot()
    for mesh in r.findall('./asset/mesh'):mesh.set('file',os.path.relpath((BASE/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for n,file in [('humerus','humerus_mounted.stl'),('ulna','forearm_mounted.stl')]:r.find(f'./asset/mesh[@name="{n}"]').set('file',file)
    ET.indent(r);tree.write(DEST/'arm.xml',encoding='utf-8',xml_declaration=True);print('mount checks',neutral,flush=True)

if __name__=='__main__':run()
