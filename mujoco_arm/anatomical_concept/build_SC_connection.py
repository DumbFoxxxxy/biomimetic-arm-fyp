"""S1 local SC ball/socket and split retention caps, with source delta evidence.
Nominal geometry only. Does not remove the bone connection to evade collision.
"""
from repair_inputs import *
from scipy.spatial import cKDTree
import os
BASE=OUT/'joint_model';DEST=BASE/'SC_connection'

def unit(x):return np.asarray(x)/np.linalg.norm(x)
def rod(a,b,r,n=24):
    a=np.asarray(a);b=np.asarray(b);z=unit(b-a);x=unit(np.cross(z,[0,0,1] if abs(z[2])<.8 else [0,1,0]));y=np.cross(z,x);return mf.Manifold.cylinder(float(np.linalg.norm(b-a)),r,r,n).transform(np.c_[np.stack([x,y,z],1),a])
def write(path,s):
    mesh=s.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3];f=np.asarray(mesh.tri_verts);tri=v[f];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.maximum(np.linalg.norm(norm,axis=1),1e-30)[:,None];dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);arr=np.zeros(len(f),dt);arr['n']=norm;arr['v']=tri;path.write_bytes(b'Independent SC connection candidate; mm; not strength validated'.ljust(80,b' ')+struct.pack('<I',len(f))+arr.tobytes());path.with_suffix('.obj').write_text('\n'.join(['# BodyParts3D-derived SC engineering candidate; CC BY-SA 2.1 Japan.']+['v '+' '.join(format(float(c),'.17g') for c in p) for p in v]+['f '+' '.join(str(int(i)+1) for i in face) for face in f])+'\n',encoding='utf8');return {'name':path.stem,'vertices':v.tolist(),'faces':f.tolist(),'volume_mm3':float(s.volume()),'components':len([c for c in s.decompose() if c.volume()>1e-6])}

def run():
    DEST.mkdir(exist_ok=True);inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));definition=json.loads((BASE/'definition.json').read_text(encoding='utf8'));S=np.array(definition['centres_source_mm']['SC'])-SHIFT;A=np.array(definition['centres_source_mm']['AC'])-SHIFT;z=unit(A-S);x=unit(np.cross(z,[0,0,1]));y=np.cross(z,x);F=np.stack([x,y,z],1)
    def local(s):return s.transform(np.c_[F,S])
    def box(lo,hi):return local(mf.Manifold.cube(np.array(hi)-lo).translate(lo))
    def load(path):
        v,f=read_obj(str(path));v,f=weld_stl(v[f]);return mf.Manifold(mf.Mesh64((v-SHIFT).astype(np.float64),f.astype(np.uint64))),v-SHIFT
    cr=next(r for r in inv['models'] if r['id']=='FJ3362');oldclav,cv=load(ROOT/cr['source']);oldman,mv=load(OUT/'repaired_supports/FJ3290.obj');ball=mf.Manifold.sphere(10,32).translate(S);pocket=mf.Manifold.sphere(22,32).translate(S)
    target=S+z*26;outer=cv[np.linalg.norm(cv-S,axis=1)>24];root=outer[np.argmin(np.linalg.norm(outer-target,axis=1))];neck=rod(S,root,4.5);clav=(oldclav-pocket)+ball+neck
    outer_shell=mf.Manifold.sphere(13,32).translate(S)-mf.Manifold.sphere(10.35,32).translate(S)
    rear=outer_shell^box([-30,-30,-30],[30,30,-.1]);ring=outer_shell^box([-30,-30,.1],[30,30,4.5]);base=oldman-mf.Manifold.sphere(18,32).translate(S)
    mount_tip=S-z*11.5;available=mv[np.linalg.norm(mv-S,axis=1)>19]
    for shift in [-5,5]:
        p=available[np.argmin(np.linalg.norm(available-(S-z*24+x*shift),axis=1))];rear=rear+rod(mount_tip+x*shift*.3,p,3)
    base=base+rear;caps=ring
    for sign in [-1,1]:
        c=S+x*sign*15;lug=rod(c+z*.1,c+z*4.5,3.7);base_lug=rod(c-z*4,c-z*.1,3.7);bridge=rod(S+x*sign*11-z*2,c-z*2,3);base=base+base_lug+bridge;caps=caps+lug
        bore=rod(c-z*14,c+z*12,1.7);base=base-bore;caps=caps-bore
    clearance=mf.Manifold.batch_boolean([caps]+[caps.translate(direction*sign*.15) for direction in [x,y,z] for sign in [-1,1]],mf.OpType.Add);base=base-clearance
    moving_clearance=mf.Manifold.batch_boolean([clav]+[clav.translate(direction*sign*.15) for direction in [x,y,z] for sign in [-1,1]],mf.OpType.Add);base=base-moving_clearance
    left=caps^box([-30,-30,-30],[-.1,30,30]);right=caps^box([.1,-30,-30],[30,30,30]);hardware={}
    for side,sign in [('left',-1),('right',1)]:
        c=S+x*sign*15;bolt=rod(c-z*10.75,c+z*5.25,1.5)+rod(c+z*5.25,c+z*8.25,2.75)
        nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[F,c-z*7.15])-rod(c-z*9,c-z*3,1.6)
        top=rod(c+z*4.65,c+z*5.15,3.5)-rod(c+z*4,c+z*6,1.7);bottom=rod(c-z*4.65,c-z*4.15,3.5)-rod(c-z*5,c-z*3,1.7)
        hardware.update({'SC_bolt_'+side:bolt,'SC_nut_'+side:nut,'SC_top_washer_'+side:top,'SC_bottom_washer_'+side:bottom})
        seat=rod(c-z*7.4,c-z*4.5,3.45)+rod(c-z*4.9,c-z*3.95,3.75);base=base-seat
    parts={};cleanup={}
    for n,s in {'clavicle_SC_ball':clav,'manubrium_socket_base':base,'SC_cap_left':left,'SC_cap_right':right,**hardware}.items():
        s=s.as_original().simplify(.01);shells=s.decompose();meaningful=[a for a in shells if abs(a.volume())>1e-10]
        if len(meaningful)==1 and len(shells)>1:cleanup[n]={'removed_shell_count':len(shells)-1,'removed_absolute_volume_mm3':sum(abs(a.volume()) for a in shells if abs(a.volume())<=1e-10)};s=meaningful[0]
        parts[n]=s
    (DEST/'mesh_cleanup.json').write_text(json.dumps({'surface_simplification_tolerance_mm':.01,'zero_volume_shell_cleanup':cleanup},indent=2),encoding='utf8');clav=parts['clavicle_SC_ball'];base=parts['manubrium_socket_base'];rows=[write(DEST/(n+'.stl'),s) for n,s in parts.items()];hits=[]
    for (an,aa),(bn,bb) in itertools.combinations(parts.items(),2):
        volume=max(0,float((aa^bb).volume()))
        if volume>.001:hits.append({'parts':[an,bn],'volume_mm3':volume})
    report={'status':'nominal SC assembly candidate; not frozen','coordinate_frame':'source raw mm minus SHIFT, same as display assets','joint_centre_source_mm':(S+SHIFT).tolist(),'ball_radius_mm':10,'socket_inner_radius_mm':10.35,'socket_outer_radius_mm':13,'neck_radius_mm':4.5,'retention_cap_depth_mm':4.5,'split_gap_mm':.2,'bolt_hole_diameter_mm':3.4,'source_overlap_mm3':float((oldclav^oldman).volume()),'neutral_part_overlaps':hits,'components':{r['name']:r['components'] for r in rows},'geometry_delta':{'clavicle_removed_mm3':float((oldclav-clav).volume()),'clavicle_added_mm3':float((clav-oldclav).volume()),'manubrium_removed_mm3':float((oldman-base).volume()),'manubrium_added_mm3':float((base-oldman).volume())},'pending':['SC full motion audit','M3 screw/nut/spacer and antirotation detail','socket support and retention load/tolerance tests','muscle origins near modified bone ends'],'physical_assembly':'insert ball into rear socket, place two split caps around neck, fasten cap tabs to rear tabs; never push ball through smaller retained aperture'};(DEST/'assembly.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');(DEST/'cad_input.json').write_text(json.dumps({'parts':rows,'units':'mm'}),encoding='utf8')
    tree=ET.parse(BASE/'arm.xml');r=tree.getroot();r.find('compiler').set('meshdir','.')
    for mesh in r.findall('./asset/mesh'):mesh.set('file',os.path.relpath((OUT/'assets'/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    for fid,name in [('FJ3362','clavicle_SC_ball'),('FJ3290','manubrium_socket_base')]:r.find(f'./asset/mesh[@name="{fid}"]').set('file',name+'.stl')
    for name in ['SC_cap_left','SC_cap_right',*hardware]:
        ET.SubElement(r.find('asset'),'mesh',name=name,file=name+'.stl',scale='.001 .001 .001');ET.SubElement(r.find('worldbody'),'geom',name=name,type='mesh',mesh=name,rgba='.2 .5 .7 1',group='3',contype='0',conaffinity='0')
    ET.indent(r);tree.write(DEST/'arm.xml',encoding='utf-8',xml_declaration=True);print('SC parts',report['components'],'overlaps',hits,flush=True)

if __name__=='__main__':run()
