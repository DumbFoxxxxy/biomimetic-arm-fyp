"""M6 exact first-sheave placement study, not a grounded transmission.

Checks actual cassette solids and a finite cable tube, rotates the real ball
terminal to its incoming tangent, and probes determinate elbow rigid geometry.
The muscle body/cap stays neutral ONLY for the placement screen: no imposed
muscle rail, contraction law, downstream support or actuation claim is made.
"""
from build_first_guide import DEST, TEMPLATE
from build_distal_swivels import ROOT, OUT, BASE, load, mf, np, intersection, unit, rod
from scipy.spatial.transform import Rotation
import json, hashlib, itertools, math
import xml.etree.ElementTree as ET

CURRENT=BASE/'swivel_terminal_distal_load'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def cable_tube(points, radius=.3, sides=12):
    p=np.asarray(points);t=np.empty_like(p);t[0]=unit(p[1]-p[0]);t[-1]=unit(p[-1]-p[-2])
    for i in range(1,len(p)-1):t[i]=unit(unit(p[i]-p[i-1])+unit(p[i+1]-p[i]))
    v=[];f=[]
    for q,d in zip(p,t):
        n=np.array([0.,0.,1.]);w=np.cross(d,n)
        v.extend(q+radius*(np.cos(a)*n+np.sin(a)*w) for a in np.arange(sides)*2*np.pi/sides)
    for i in range(len(p)-1):
        for j in range(sides):
            a=i*sides+j;b=i*sides+(j+1)%sides;c=b+sides;d=a+sides;f.extend([[a,b,c],[a,c,d]])
    for j in range(1,sides-1):f.extend([[0,j+1,j],[(len(p)-1)*sides,(len(p)-1)*sides+j,(len(p)-1)*sides+j+1]])
    s=mf.Manifold(mf.Mesh64(np.asarray(v),np.asarray(f,dtype=np.uint64)))
    assert s.status()==mf.Error.NoError and s.volume()>0
    return s

def route(L,beta,gamma):
    """Coordinates about sheave origin; actual circle tangencies, clockwise wrap."""
    r=8.005;b=np.radians(beta);g=np.radians(gamma);p=-L*np.array([np.sin(b),np.cos(b)])
    tangent=(r*r/(L*L))*p-r*np.sqrt(L*L-r*r)/(L*L)*np.array([-p[1],p[0]])
    incoming=unit(tangent-p);start=p+incoming*3.5
    a=math.atan2(tangent[1],tangent[0]);wrap=(a-g)%(2*np.pi)
    angles=np.linspace(a,a-wrap,int(np.ceil(np.degrees(wrap)))+1)
    arc=np.c_[r*np.cos(angles),r*np.sin(angles)]
    outgoing=np.array([np.sin(g),-np.cos(g)])
    points=np.vstack([start,arc,arc[-1]+outgoing*40])
    return cable_tube(np.c_[points,np.zeros(len(points))]),np.r_[incoming,0],np.r_[p,0],np.degrees(wrap)

def hits(s,world):
    return [{'part':k,'volume_mm3':v} for k,b in world.items() if (v:=intersection(s,b))>.001]

def run():
    design=json.loads((TEMPLATE/'design.json').read_text());parts={k:load(TEMPLATE/(k+'.obj'))[0] for k in design['parts']}
    source_hashes={str(TEMPLATE/(k+'.obj')):digest(TEMPLATE/(k+'.obj')) for k in parts}
    local={'scope':'cassette solid geometry and nominal 0.6mm cable only; no material, strength, wear or cable minimum bend rating', 'source_sha256':source_hashes,'neutral_interferences':design['neutral_interferences']}
    # One full wheel revolution, all other actual parts retained in the check.
    obstacles={k:s for k,s in parts.items() if k!='grooved_wheel'}
    rotation_hits=[]
    for angle in range(361):
        R=Rotation.from_rotvec([0,0,np.radians(angle)]).as_matrix()
        h=hits(parts['grooved_wheel'].transform(np.c_[R,np.zeros(3)]),obstacles)
        if h:rotation_hits.append({'angle_deg':angle,'hits':h})
    local['wheel_rotation']={'step_deg':1,'forward_samples':361,'return_same_configurations':True,'hits':rotation_hits}
    testwire,_,_,wrap=route(80,20,45)
    local['nominal_wrap_example']={'source_distance_mm':80,'offset_angle_deg':20,'outgoing_angle_deg':45,'wrap_deg':wrap,'arc_step_max_deg':1,'hits':hits(testwire,parts)}
    local['clamp_load_path']='axle head -> front washer -> fork cheek -> sleeve -> rear cheek -> washer -> nut; nominal geometry, preload/creep unverified'
    (DEST/'cassette_checks.json').write_text(json.dumps(local,indent=2),encoding='utf8')
    assert not rotation_hits and not local['nominal_wrap_example']['hits'] and not design['neutral_interferences']
    print('cassette: wheel 361 samples and real tangent/arc cable clear',flush=True)

    a=json.loads((CURRENT/'assembly.json').read_text());world={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in a['instances']}
    kinds={r['key']:r['kind'] for r in a['instances']}
    for r in a['instances']:source_hashes[str(ROOT/r['file'])]=digest(ROOT/r['file'])
    support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in support['models']:
        if r['kind']!='bone_reference':continue
        if r['id'] in world:continue
        # Use the already documented, closed computational reference copies;
        # never silently cap the original BodyParts3D meshes during this screen.
        p=OUT/'repaired_supports'/(r['id']+'.obj');world['support_'+r['id']]=load(p,raw=True)[0];kinds['support_'+r['id']]='reference_bone';source_hashes[str(p)]=digest(p)
    definition=json.loads((OUT/'joint_model/definition.json').read_text());E=np.array(definition['centres_source_mm']['elbow'])-np.array(definition['source_coordinate_shift_mm'])
    xml=OUT/'elbow_muscles/mounted_elbow_candidate/arm.xml';hinge=ET.parse(xml).find(".//joint[@name='elbow_hinge']");axis=np.fromstring(hinge.get('axis'),sep=' ')
    # Loose old tendon stops are packing references, included in neutral checks
    # but are not treated as bolted moving parts in the sweep.
    moving={r['key']:world[r['key']] for r in a['instances'] if r['owner']=='elbow' and not r['key'].endswith('_tendon_stop')}
    cassette=mf.Manifold.batch_boolean(list(parts.values()),mf.OpType.Add)
    results=[]
    for fid in ['FJ1486','FJ1477']:
        d=json.loads((CURRENT/fid/'design.json').read_text());C=np.array(d['selected']['ball_centre_display_mm']);u=unit(np.array(d['selected']['neutral_terminal_axis']))
        v=unit(axis-u*(axis@u));w=np.cross(u,v);ballkey=fid+'_swivel_ball_crimp_blank';ball=world[ballkey];other={k:s for k,s in world.items() if k!=ballkey}
        trials=[];feasible=[]
        for L,beta,azimuth in itertools.product([30,50,80,110,150],[0,20],range(0,360,45)):
            wire,inc,p,wrap=route(L,beta,45);phi=np.radians(azimuth);x=v*np.cos(phi)+w*np.sin(phi);F=np.stack([x,u,np.cross(x,u)],axis=1);Q=C-F@p;T=np.c_[F,Q]
            cassette_world=cassette.transform(T);wire_world=wire.transform(T);direction=F@inc;tilt=np.degrees(np.arccos(np.clip(direction@u,-1,1)))
            row={'source_distance_mm':L,'offset_angle_deg':beta,'azimuth_deg':azimuth,'outgoing_angle_deg':45,'terminal_tilt_deg':float(tilt),'wrap_deg':wrap,'centre_display_mm':Q.tolist(),'frame_columns':F.tolist()}
            if tilt>25:row['rejected']='terminal_angle';trials.append(row);continue
            rotvec=np.cross(u,direction);norm=np.linalg.norm(rotvec);R=Rotation.from_rotvec(unit(rotvec)*np.radians(tilt)).as_matrix() if norm>1e-10 else np.eye(3)
            rotated_ball=ball.transform(np.c_[R,C-R@C]);checks={}
            for label,s,against in [('cassette',cassette_world,world),('cable',wire_world,{**other,ballkey:rotated_ball}),('articulated_terminal',rotated_ball,other)]:
                h=hits(s,against)
                if h:checks[label]=h
            # The wire must also clear the real guide; no unsupported kink.
            h=hits(wire_world,{'guide':cassette_world})
            if h:checks['cable_guide']=h
            row['neutral_interferences']=checks;row['neutral_clear']=not checks;trials.append(row)
            if not checks:feasible.append((row,cassette_world,wire_world,rotated_ball))
        print(fid,'neutral candidates',len(trials),'clear',len(feasible),flush=True)
        # Sweep every neutral-clear placement in this finite search family.
        sweeps=[]
        for row,s,wire,terminal in feasible:
            first=None;last=0
            for angle in range(151):
                R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E]
                h=hits(s,{k:b.transform(T) for k,b in moving.items()})
                if h:
                    lo=max(0,angle-1);hi=float(angle)
                    while hi-lo>.10000001:
                        mid=(lo+hi)/2;R=Rotation.from_rotvec(axis*np.radians(mid)).as_matrix()
                        if hits(s,{k:b.transform(np.c_[R,E-R@E]) for k,b in moving.items()}):hi=mid
                        else:lo=mid
                    last=lo;first={'upper_angle_deg':hi,'hits_at_integer_angle':h};break
                last=angle
            sweeps.append({**row,'continuous_forward_prefix_deg':last,'first_rigid_limit':first,'angular_step_deg':1,'boundary_interval_max_deg':.1,'return_uses_identical_prescribed_poses':True})
        selected=max(sweeps,key=lambda r:(r['continuous_forward_prefix_deg'],-r['source_distance_mm'])) if sweeps else None
        results.append({'muscle':fid,'trials':trials,'neutral_clear_count':len(feasible),'sweeps':sweeps,'best_placement_only':selected})
        print(fid,'best guide-only elbow prefix',selected['continuous_forward_prefix_deg'] if selected else None,flush=True)
    report={'scope':'M6 UNGROUNDED first-sheave placement screen. All manufactured neutral bodies and bony references checked; finite 40mm outgoing corridor has no downstream constraint. Only determinate bolted elbow rigid bodies move in sweep; no cap/body motion, muscle contraction or assembled transmission proof.','modules':results,'source_sha256':source_hashes,'moving_rigid_keys':list(moving),'mount_frame_designed':False,'full_transmission_exists':False,'frozen':False}
    (DEST/'first_guide_placement_screen.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    assert all(digest(__import__('pathlib').Path(p))==h for p,h in source_hashes.items())

if __name__=='__main__':run()
