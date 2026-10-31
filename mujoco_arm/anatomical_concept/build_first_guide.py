"""M6 physical turnaround cassette, local mm coordinates, not yet grounded.

LocalY is the wrap side, localZ is the axle; the keyed mount is behind the wrap.
The shoulder sleeve carries nominal axle clamp load without pinching the wheel.
"""
from build_distal_swivels import ROOT,OUT,load,rod,unit,mf,np,write,intersection
import json,itertools
DEST=OUT/'elbow_guides'
TEMPLATE=DEST/'template'

def annular_profile(profile,n=48):
    vertices=np.array([[r*np.cos(a),r*np.sin(a),z] for r,z in profile for a in np.arange(n)*2*np.pi/n]);faces=[]
    for ring in range(len(profile)):
        nxt=(ring+1)%len(profile)
        for j in range(n):
            a=ring*n+j;b=ring*n+(j+1)%n;c=nxt*n+(j+1)%n;d=nxt*n+j;faces += [[a,b,c],[a,c,d]]
    s=mf.Manifold(mf.Mesh64(vertices,np.array(faces,dtype=np.uint64)));assert s.status()==mf.Error.NoError and s.volume()>0;return s

def axle_bolt():
    n=48;v=[];f=[]
    for r,z in [(1.5,-10.6),(1.5,5.4),(2.75,5.4),(2.75,8.4)]:v += [[r*np.cos(a),r*np.sin(a),z] for a in np.arange(n)*2*np.pi/n]
    v += [[0,0,-10.6],[0,0,8.4]]
    for ring in range(3):
        for j in range(n):
            a=ring*n+j;b=ring*n+(j+1)%n;c=(ring+1)*n+(j+1)%n;d=(ring+1)*n+j;f += [[a,b,c],[a,c,d]]
    for j in range(n):f += [[(j+1)%n,j,4*n],[3*n+j,3*n+(j+1)%n,4*n+1]]
    return mf.Manifold(mf.Mesh64(np.array(v),np.array(f,dtype=np.uint64)))

def run():
    TEMPLATE.mkdir(parents=True,exist_ok=True)
    wheel=annular_profile([(9.5,-2.25),(9.5,-.6),(7.7,-.6),(7.7,.6),(9.5,.6),(9.5,2.25),(2.65,2.25),(2.65,-2.25)])
    sleeve=annular_profile([(2.5,-2.4),(2.5,2.4),(1.6,2.4),(1.6,-2.4)])
    fork=rod([0,0,2.4],[0,0,4.9],4,48)+rod([0,0,-4.9],[0,0,-2.4],4,48)
    fork=fork+mf.Manifold.cube([4.4,12,2.5]).translate([-2.2,-13,2.4])+mf.Manifold.cube([4.4,12,2.5]).translate([-2.2,-13,-4.9])
    fork=fork+rod([0,-13,-4.9],[0,-13,4.9],3,48)+mf.Manifold.cube([6,6,10]).translate([-3,-16,-5])
    fork=(fork-rod([0,0,-12],[0,0,12],1.7,48)-rod([0,-13,-12],[0,-13,12],1.7,48)).as_original().simplify(.001)
    washer_front=annular_profile([(3.5,4.9),(3.5,5.4),(1.6,5.4),(1.6,4.9)])
    washer_back=annular_profile([(3.5,-5.4),(3.5,-4.9),(1.6,-4.9),(1.6,-5.4)])
    nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,0,-7.8])-rod([0,0,-9],[0,0,-4],1.6,48)
    parts={'cassette_fork':fork,'grooved_wheel':wheel,'shoulder_sleeve':sleeve,'M3x16_axle':axle_bolt(),'front_washer':washer_front,'rear_washer':washer_back,'M3_nut':nut};rows=[]
    for name,s in parts.items():
        assert s.status()==mf.Error.NoError;row=write(TEMPLATE/(name+'.stl'),s);assert row['components']==1;rows.append(row)
    overlaps=[]
    for (a,sa),(b,sb) in itertools.combinations(parts.items(),2):
        volume=intersection(sa,sb)
        if volume>.001:overlaps.append({'pair':[a,b],'volume_mm3':volume})
    design={'scope':'M6 physical first-turn cassette geometry only, requires actual keyed ground frame','units':'mm','parts':{r['name']:r['components'] for r in rows},'wire_radius_mm':.3,'wire_trace_radius_mm':8.005,'groove_root_radius_mm':7.7,'groove_half_width_mm':.6,'wheel_flange_radius_mm':9.5,'wheel_half_width_mm':2.25,'wheel_bore_radius_mm':2.65,'sleeve_outer_radius_mm':2.5,'sleeve_inner_radius_mm':1.6,'sleeve_end_faces_mm':[-2.4,2.4],'nominal_radial_bearing_clearance_mm':.15,'nominal_axial_bearing_clearance_each_end_mm':.15,'keyed_mount_local_centre_mm':[0,-13,0],'keyed_mount_size_mm':[6,6,10],'mount_bore_radius_mm':1.7,'mount_frame_and_hardware_designed':False,'metal_hardware':'nominal envelopes; actual fastener/bushing material and fit unselected','neutral_interferences':overlaps,'bearing_load_or_efficiency_verified':False,'frozen':False}
    (TEMPLATE/'design.json').write_text(json.dumps(design,indent=2),encoding='utf8');(TEMPLATE/'cad_input.json').write_text(json.dumps({'units':'mm','parts':rows}),encoding='utf8');print('M6 guide pieces',len(parts),'overlaps',overlaps,flush=True)
if __name__=='__main__':run()
