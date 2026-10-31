"""M8 physical two-axis receiving crimp terminals on existing lever eyelets.

Yaw about an elbow-parallel shoulder sleeve plus pitch about a transverse
shoulder sleeve; neither motion is replaced by a fictitious free cable kink.
No upstream tendon route or crimp/preload/strength claim is made here.
"""
from screen_first_guides import ROOT,OUT,BASE,load,mf,np,rod,unit,intersection
from build_first_guide import annular_profile,write
from scipy.spatial.transform import Rotation
import json,hashlib,itertools

PARENT=OUT/'elbow_guides/ground_frame_swept'
DEST=OUT/'forearm_receivers'

def bolt(length,head_z,shaft_r,head_r,head_height,n=32):
    vertices=[];faces=[]
    for r,z in [(shaft_r,head_z-length),(shaft_r,head_z),(head_r,head_z),(head_r,head_z+head_height)]:vertices.extend([r*np.cos(a),r*np.sin(a),z] for a in np.arange(n)*2*np.pi/n)
    vertices.extend([[0,0,head_z-length],[0,0,head_z+head_height]])
    for i in range(3):
        for j in range(n):
            a=i*n+j;b=i*n+(j+1)%n;c=b+n;d=a+n;faces.extend([[a,b,c],[a,c,d]])
    for j in range(n):faces.extend([[(j+1)%n,j,4*n],[3*n+j,3*n+(j+1)%n,4*n+1]])
    s=mf.Manifold(mf.Mesh64(np.array(vertices),np.array(faces,dtype=np.uint64)));assert s.status()==mf.Error.NoError;return s

def construct(standoff=0):
    h=5.5+standoff
    sleeve=annular_profile([(3.8,2.3),(3.8,3.8+standoff),(2.5,3.8+standoff),(2.5,7.2+standoff),(1.6,7.2+standoff),(1.6,2.3)],32)
    rotor=annular_profile([(4.5,4+standoff),(4.5,7+standoff),(2.65,7+standoff),(2.65,4+standoff)],32)
    for sign in [-1,1]:
        a,b=sorted([sign*2.4,sign*4.6]);ear=rod([8,a,h],[8,b,h],3.4,32);arm=rod([2.5,sign*2,h],[8,sign*3.5,h],1.6,24);rotor=rotor+ear+arm
    rotor=rotor^mf.Manifold.cube([60,9.2,60]).translate([-20,-4.6,-15])
    rotor=(rotor-rod([8,-2.4,h],[8,2.4,h],2.65,32)-rod([8,-8,h],[8,8,h],1.1,32)-rod([0,0,-3],[0,0,standoff+12],2.65,32)).as_original().simplify(.001)
    hub=(rod([8,-2,h],[8,2,h],2.5,32)+rod([9.5,0,h],[12.2,0,h],1.2,32)+rod([12,0,h],[17.8,0,h],1.6,32)-rod([8,-4,h],[8,4,h],1.55,32)-rod([12.1,0,h],[19,0,h],.45,32)).as_original().simplify(.001)
    # Local pitch hardware is built along Z then mapped onto +Y at X=8.
    S=np.array([[1.,0,0],[0,0,1],[0,-1,0]]);T=np.c_[S,[8,0,h]]
    parts={'yaw_sleeve':sleeve,'rotor_yoke':rotor,'pitch_crimp_blank':hub,
      'pitch_sleeve':annular_profile([(1.4,-2.4),(1.4,2.4),(1.1,2.4),(1.1,-2.4)],32).transform(T),
      'pitch_M2x12':bolt(12,5.,1.,1.9,1.6).transform(T),
      'pitch_front_washer':annular_profile([(2.5,4.6),(2.5,5),(1.1,5),(1.1,4.6)],32).transform(T),
      'pitch_back_washer':annular_profile([(2.5,-5),(2.5,-4.6),(1.1,-4.6),(1.1,-5)],32).transform(T),
      'pitch_nut':(mf.Manifold.cylinder(1.6,4/np.sqrt(3),4/np.sqrt(3),6).translate([0,0,-6.6])-rod([0,0,-8],[0,0,-4],1.05,32)).transform(T),
      'yaw_M3':bolt(next(L for L in [16,20,25,30,35] if L>=13.4+standoff),7.7+standoff,1.5,2.75,3),
      'yaw_front_washer':annular_profile([(3.5,7.2+standoff),(3.5,7.7+standoff),(1.6,7.7+standoff),(1.6,7.2+standoff)],32),
      'yaw_base_washer':annular_profile([(3.5,2),(3.5,2.3),(1.6,2.3),(1.6,2)],32),
      'yaw_back_washer':annular_profile([(3.5,-2.3),(3.5,-2),(1.6,-2),(1.6,-2.3)],32),
      'yaw_nut':mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,0,-4.7])-rod([0,0,-6],[0,0,0],1.6,32)}
    return parts

YAW_PARTS={'rotor_yoke','pitch_sleeve','pitch_M2x12','pitch_front_washer','pitch_back_washer','pitch_nut','pitch_crimp_blank'}

def articulate(parts,yaw_deg,pitch_deg,standoff):
    R=Rotation.from_rotvec([0,0,np.radians(yaw_deg)]).as_matrix();P=Rotation.from_rotvec([0,-np.radians(pitch_deg),0]).as_matrix();H=np.array([8,0,5.5+standoff]);out={}
    for key,s in parts.items():
        if key=='pitch_crimp_blank':s=s.transform(np.c_[P,H-P@H])
        if key in YAW_PARTS:s=s.transform(np.c_[R,np.zeros(3)])
        out[key]=s
    return out

def run():
    DEST.mkdir(exist_ok=True);reports=[]
    for fid,standoff in [('FJ1486',6),('FJ1477',11)]:
        folder=DEST/fid;folder.mkdir(exist_ok=True);parts=construct(standoff);rows=[]
        for key,s in parts.items():
            assert s.status()==mf.Error.NoError;row=write(folder/(key+'.stl'),s);assert row['components']==1;rows.append(row)
        overlaps=[]
        for (a,sa),(b,sb) in itertools.combinations(parts.items(),2):
            v=intersection(sa,sb)
            if v>.001:overlaps.append({'pair':[a,b],'volume_mm3':v})
        d={'scope':'M8 real 2-axis receiving-terminal manufacturing geometry, local mm; no installed upstream route or strength/crimp proof','muscle':fid,'parts':[r['name'] for r in rows],'standoff_mm':standoff,'yaw_axis':[0,0,1],'yaw_centre':[0,0,5.5+standoff],'pitch_axis_at_zero_yaw':[0,1,0],'pitch_centre_at_zero_yaw':[8,0,5.5+standoff],'wire_outlet_at_zero_pose':[17.8,0,5.5+standoff],'wire_bore_radius_mm':.45,'wire_radius_mm':.3,'yaw_radial_clearance_mm':.15,'yaw_axial_clearance_each_mm':.2,'pitch_radial_clearance_mm':.15,'pitch_axial_clearance_each_mm':.4,'neutral_interferences':overlaps,'force_path':'wire crimp blank -> pitch shoulder sleeve/M2 pin -> rotor yoke -> yaw shoulder sleeve/M3 pin -> lever eyelet -> ulna mounting bolts','not_verified':['crimp grip','nominal fastener/thread/supplier selection','materials','preload/friction/wear','full upstream routing','muscle deformation and arm drive'],'frozen':False}
        (folder/'design.json').write_text(json.dumps(d,indent=2),encoding='utf8');(folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':rows}),encoding='utf8');reports.append(d);print(fid,'pieces',len(parts),'neutral hits',overlaps,flush=True)
    (DEST/'modules.json').write_text(json.dumps({'modules':[d['muscle'] for d in reports],'scope':'two receiving mechanism candidates only','frozen':False},indent=2),encoding='utf8')

if __name__=='__main__':run()
