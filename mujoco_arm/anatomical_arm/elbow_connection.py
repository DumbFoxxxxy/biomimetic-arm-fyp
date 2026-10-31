"""Outboard elbow bearing bridge, preserving joint surfaces except axle bore.
Independent engineering candidate; does not overwrite anatomical A1.
"""
from build import *
from scipy.spatial.transform import Rotation
DEST=OUT/'elbow_connection'

def load(n):
    v,f=read_stl(str(OUT/'assets'/(n+'.stl')))
    return mf.Manifold(mf.Mesh(v.astype(np.float32),f.astype(np.uint32)))

def nearest(s,p):
    v=np.asarray(s.to_mesh64().vert_properties)[:,:3];return v[np.argmin(np.linalg.norm(v-p,axis=1))]

def build_candidate(offset):
    a=json.loads((OUT/'source_audit.json').read_text(encoding='utf8'));lm=a['landmarks_mm'];E=np.array(lm['E']);axis=unit(lm['elbow_axis']);z=unit(E-lm['GH']);h=load('humerus');u=load('ulna');r=load('radius')
    bore=rod(E-axis*60,E+axis*60,3.15,32);hum=h-bore
    bridges=[]
    for sign in [-1,1]:
        c=E+axis*offset*sign;ring=rod(c-axis*4,c+axis*4,8,32)-rod(c-axis*5,c+axis*5,3.3,32)
        root=nearest(u,E+z*65+axis*sign*8)
        # Keep two beams outboard; enter the annular material radially.
        anchor=c+z*6;corner=root+axis*sign*(offset-8)
        bridges.extend([ring,rod(anchor,corner,3.5,16),rod(corner,root,3.5,16)])
    p=nearest(u,E+z*130);q=nearest(r,E+z*130)
    fore=union([u,r,rod(p,q,4,16),*bridges])-bore
    pin=rod(E-axis*(offset+10),E+axis*(offset+10),3,32)
    head=rod(E-axis*(offset+13),E-axis*(offset+10),6,32)
    # Positive retention: cross-drilled end for a cotter pin, nominal 1.6 mm.
    end=E+axis*(offset+8);cotter_axis=unit(np.cross(axis,z));pin=(pin+head)-rod(end-cotter_axis*5,end+cotter_axis*5,.8,16)
    cotter=union([rod(end-cotter_axis*5,end+cotter_axis*5,.7,16),ball(end-cotter_axis*5,1.2),rod(end+cotter_axis*5,end+cotter_axis*5+axis*3,.7,16)])
    parts={'humerus_bored':hum,'forearm_bridged':fore,'retained_axle':pin,'cotter_pin':cotter}
    for sign in [-1,1]:
        c=E+axis*sign*(offset+4.6);parts['thrust_washer_'+str(sign)]=rod(c-axis*.5,c+axis*.5,7,32)-rod(c-axis,c+axis,3.2,32)
        parts['forearm_bridged']=parts['forearm_bridged']-rod(c-axis*.65,c+axis*.65,7.2,32)
    return parts,E,axis

def audit(parts,E,axis):
    h=parts['humerus_bored'];f=parts['forearm_bridged'];pin=parts['retained_axle'];samples=[];first=None
    def accessory_hits(R):
        T=np.c_[R,E-R@E];hits=[]
        for mn in ['forearm_bridged','thrust_washer_-1','thrust_washer_1']:
            moving=parts[mn].transform(T)
            for sn in ['humerus_bored','retained_axle','cotter_pin']:
                vol=float((moving^parts[sn]).volume())
                if vol>.001:hits.append({'parts':[mn,sn],'volume_mm3':vol})
        return hits
    for deg in range(146):
        R=Rotation.from_rotvec(axis*np.radians(deg)).as_matrix();m=f.transform(np.c_[R,E-R@E]);vol=max(0,float((h^m).volume()));pv=max(0,float((pin^m).volume()));accessories=accessory_hits(R);hit=bool(accessories)
        samples.append({'angle_deg':deg,'humerus_forearm_overlap_mm3':vol,'axle_forearm_overlap_mm3':pv,'all_moving_fixed_component_hits':accessories,'minimum_gap_mm_capped_2':float(h.min_gap(m,2))})
        if hit and first is None:first=deg
    safe=145 if first is None else max(0,first-1);boundary=None
    if first and first>0:
        lo=first-1;hi=first
        while hi-lo>.1:
            mid=(lo+hi)/2;R=Rotation.from_rotvec(axis*np.radians(mid)).as_matrix();m=f.transform(np.c_[R,E-R@E])
            if accessory_hits(R):hi=mid
            else:lo=mid
        safe=lo;boundary=hi
    reverse=[]
    for deg in reversed(range(int(safe)+1)):
        R=Rotation.from_rotvec(axis*np.radians(deg)).as_matrix();m=f.transform(np.c_[R,E-R@E]);reverse.append({'angle_deg':deg,'humerus_forearm_overlap_mm3':max(0,float((h^m).volume())),'axle_forearm_overlap_mm3':max(0,float((pin^m).volume())),'all_moving_fixed_component_hits':accessory_hits(R)})
    neutral_pairs=[]
    for (an,aa),(bn,bb) in itertools.combinations(parts.items(),2):
        vol=float((aa^bb).volume())
        if vol>.001:neutral_pairs.append({'parts':[an,bn],'overlap_mm3':vol})
    return {'continuous_prefix_deg':safe,'first_collision_deg':boundary if boundary else first,'forward_samples':samples,'return_samples':reverse,'neutral_component_overlaps':neutral_pairs,'neutral_axle_humerus_overlap_mm3':float((h^pin).volume()),'components':{n:len([x for x in s.decompose() if x.volume()>1e-6]) for n,s in parts.items()}}

def run():
    DEST.mkdir(exist_ok=True);candidates=[];best=None
    for offset in [30.,36.,42.]:
        parts,E,axis=build_candidate(offset);result=audit(parts,E,axis);result['bearing_offset_each_side_mm']=offset;candidates.append(result)
        print('elbow bridge',offset,'ROM',result['continuous_prefix_deg'],'components',result['components'],flush=True)
        score=(all(v==1 for v in result['components'].values()),result['continuous_prefix_deg'],-offset)
        if best is None or score>best[0]:best=(score,parts,result)
        if score[0] and result['continuous_prefix_deg']==145:break
    _,parts,result=best;cad=[]
    for n,s in parts.items():cad.append({'name':n,**stl(DEST/(n+'.stl'),s)})
    dump(DEST/'cad_input.json',{'parts':cad,'units':'mm','format':'elbow manufacturing geometry candidate'})
    dump(DEST/'audit.json',{'selected':result,'candidates':candidates,'retention':'integral axle head; 1.4mm cotter in 1.6mm cross-hole, enlarged head and tail bent AFTER insertion; separate 1mm thrust washers','bore_diameter_mm':6.3,'bearing_bore_diameter_mm':6.6,'axle_diameter_mm':6,'validation':'geometry only; no strength/tolerance/wear acceptance; attachment to full arm still requires audit. Axle lock to humerus and washer axial play need mechanical review.','source_anatomical_A1_unchanged':True})

if __name__=='__main__':run()
