"""Local sternum-mounted two-port reaction frame; no engineering backboard."""
from itertools import combinations
import hashlib
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,export,rod,mesh_geom
from build_channel_first_inventory import provenance


def line_crossings(tri,x,z):
    origin=np.array([x,-260.,z]);direction=np.array([0.,1.,0.])
    e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];h=np.cross(np.broadcast_to(direction,e2.shape),e2)
    determinant=np.sum(e1*h,axis=1);valid=abs(determinant)>1e-10
    inverse=np.divide(1.,determinant,out=np.zeros_like(determinant),where=valid)
    s=origin-tri[:,0];u=inverse*np.sum(s*h,axis=1);q=np.cross(s,e1);v=inverse*(q@direction)
    distance=inverse*np.sum(e2*q,axis=1)
    valid&=(u>=-1e-8)&(v>=-1e-8)&(u+v<=1+1e-8)
    hits=np.unique(np.round(distance[valid],8));assert len(hits)==2,hits
    return origin+direction*hits[0],origin+direction*hits[1]


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'elbow_supply_interface/assembly.json';parent=read(parent_file)
    dest=base/'thorax_reaction_mount';dest.mkdir(exist_ok=True)
    solids={k:load_component(v) for k,v in parent['sources'].items()};sternum_key='thorax_reference:FJ3178'
    old=solids[sternum_key];vertices,faces=mesh_geom.manifold_mesh_to_arrays(old.to_mesh64())
    # Retain a tiny disconnected source fragment by adding a bounded bridge;
    # never delete an input island to satisfy a connected-material check.
    native=sorted(old.decompose(),key=lambda s:s.volume(),reverse=True)
    connected=old;source_repairs=[]
    if len(native)>1:
        mv,mf0=mesh_geom.manifold_mesh_to_arrays(native[0].to_mesh64());tri=mv[mf0]
        for fragment in native[1:]:
            bounds=np.array(fragment.bounding_box());centre=(bounds[:3]+bounds[3:])/2
            distance,closest=mesh_geom._closest_all(centre,tri);j=int(np.argmin(distance))
            signed,normal=mesh_geom.signed_distance_and_normal(centre,tri)
            assert fragment.volume()<.01 and abs(signed)<.5,'Source repair exceeds local microfragment policy'
            bridge=rod(centre,closest[j]-normal*.4,.3);untrimmed=bridge;reliefs=[]
            for key,obstacle in solids.items():
                if key==sternum_key:continue
                volume=float(A.I.intersection_volume(bridge,obstacle))
                if volume>.001:
                    bridge-=obstacle;reliefs.append(dict(source_key=key,removed_overlap_mm3=volume))
            connected+=bridge
            source_repairs.append(dict(fragment_volume_mm3=float(fragment.volume()),gap_mm=abs(float(signed)),
                fragment_centre_source_mm=centre.tolist(),surface_point_source_mm=closest[j].tolist(),
                bridge_radius_mm=.3,bridge_new_material_mm3=float((bridge-old).volume()),
                neighbour_reliefs=reliefs))
    mf=A.I._manifold();cube=mf.Manifold.cube
    box=lambda lo,hi:cube(np.array(hi)-lo).translate(lo)
    normal=np.array([0.,-1.,0.]);mounts=[];stock=[];holes=[];parts={}
    frame=box([-10.,-222.,1220.],[10.,-218.,1282.])
    for i,z in enumerate((1230.,1260.)):
        front,back=line_crossings(vertices[faces],0.,z)
        seat_front=front+normal*2.;seat_back=back-normal*2.
        stock.extend([rod(front-normal*3.,seat_front,5.),rod(back+normal*3.,seat_back,5.)])
        hole=rod([0.,-260.,z],[0.,-140.,z],1.7);holes.append(hole)
        parts[f'thorax_mount:spacer_{i}']=rod([0.,-218.,z],seat_front,5.)-hole
        head=np.array([0.,-222.5,z]);length=float(5*np.ceil(((head-seat_back)@normal+3.9)/5.))
        parts[f'thorax_mount:bolt_{i}']=rod(head,head-normal*length,1.5)+rod(head,head+normal*3.,2.75)
        parts[f'thorax_mount:front_washer_{i}']=rod([0.,-222.,z],head,3.5)-hole
        parts[f'thorax_mount:rear_washer_{i}']=rod(seat_back,seat_back-normal*.5,3.5)-hole
        H=np.stack([[1.,0.,0.],[0.,0.,-1.],[0.,1.,0.]],axis=1)
        nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[H,seat_back-normal*.5])
        parts[f'thorax_mount:nut_{i}']=nut-hole
        mounts.append(dict(z_source_mm=z,front_source_mm=front.tolist(),back_source_mm=back.tolist(),
            anterior_seat_source_mm=seat_front.tolist(),posterior_seat_source_mm=seat_back.tolist(),
            bolt_axis_source=normal.tolist(),nominal_M3_length_mm=length,boss_radius_mm=5.,bore_radius_mm=1.7,
            source_thickness_along_axis_mm=float((front-back)@normal)))
    sternum=(connected+mf.Manifold.batch_boolean(stock,mf.OpType.Add))-mf.Manifold.batch_boolean(holes,mf.OpType.Add)
    frame-=mf.Manifold.batch_boolean(holes,mf.OpType.Add)
    ports=[];stub_solids={}
    for head,z in [('FJ1486',1236.),('FJ1477',1272.)]:
        C=np.array([-50.,-235.,z]);axis=np.array([-1.,0.,0.]);F=np.stack([[0.,0.,1.],[0.,1.,0.],axis],axis=1)
        assert np.linalg.det(F)>.999999
        pt=lambda a,b,t:C+F@np.array([a,b,t])
        cyl=lambda a,b,t0,t1,r:rod(pt(a,b,t0),pt(a,b,t1),r)
        body=cyl(0,0,-.5,11.2,4.)
        cuts=[cyl(0,0,-.6,4.,.85),cyl(0,0,4.,10.,2.1),cyl(0,0,10.,11.3,3.1)]
        plate=cyl(0,0,11.2,13.2,4.)-cyl(0,0,11.1,13.3,1.7)
        ear_positions=[]
        for j,sign in enumerate((-1.,1.)):
            a,b=sign*9/np.sqrt(2),-9/np.sqrt(2);radial=np.array([a,b])/9.;tangent=np.array([radial[1],-radial[0]])
            G=F@np.array([[tangent[0],radial[0],0],[tangent[1],radial[1],0],[0,0,1]])
            def earbox(lo,hi):return cube(np.array(hi)-lo).translate(lo).transform(np.c_[G,pt(a,b,0)])
            body+=cyl(a,b,-.5,11.2,3.2)+rod(pt(*(radial*3),5.85),pt(a,b,5.85),1.6)
            cuts.append(earbox([-2.1,-2.6,-.25],[2.1,5.,1.75]))
            hole=cyl(a,b,-.6,16.,1.15);cuts.append(hole)
            plate+=rod(pt(*(radial*3),12.2),pt(a,b,12.2),1.)+cyl(a,b,11.2,13.2,3.2)
            plate-=hole
            stem='thorax_port:'+head
            parts[stem+f'_bolt_{j}']=cyl(a,b,-.5,13.5,1.)+cyl(a,b,13.5,15.5,1.9)
            parts[stem+f'_washer_{j}']=cyl(a,b,13.2,13.5,2.5)-hole
            nut=mf.Manifold.cylinder(1.6,4./np.sqrt(3),4./np.sqrt(3),6).rotate([0,0,30]).transform(np.c_[G,pt(a,b,.15)])
            parts[stem+f'_nut_{j}']=nut-cyl(a,b,.05,1.85,1.05)
            ear_positions.append([a,b])
        body-=mf.Manifold.batch_boolean(cuts,mf.OpType.Add)
        # Open bracket below the wire corridor, not a plate behind the arm.
        frame+=body+rod([-8.,-220.,z-8.],C+[0.,0.,-8.],2.5)+rod(C+[0.,0.,-8.],C+[0.,0.,-3.],2.)
        ferrule=(cyl(0,0,4.,10.,2.)+cyl(0,0,10.,11.,3.))-cyl(0,0,3.9,11.1,.65)-cyl(0,0,7.,11.1,1.55)
        parts['thorax_port:'+head+'_ferrule']=ferrule
        parts['thorax_port:'+head+'_retainer']=plate
        stub_solids[head]=cyl(0,0,7.,24.,1.5)-cyl(0,0,6.9,24.1,.6)
        ports.append(dict(id=head,owner='thorax',frame_columns=F.tolist(),frame_origin_source_mm=C.tolist(),
            ferrule_bearing_source_mm=pt(0,0,4.).tolist(),sheath_reaction_seat_source_mm=pt(0,0,7.).tolist(),
            ferrule_outlet_source_mm=pt(0,0,11.).tolist(),outward_axis_source=axis.tolist(),
            reference_only_stub_tip_source_mm=pt(0,0,24.).tolist(),reference_tip_fixed=False,
            wire_inboard_probe_source_mm=pt(0,0,-20.).tolist(),ear_centres_local_mm=ear_positions,
            reaction_chain='ferrule inner step -> ferrule bottom -> integral port/frame -> M3 front rail/spacers -> modified sternum bearing bosses -> sternum; rib/spine/base load path pending',
            drive_end='No motor/spool selected or mounted; inboard wire corridor is a probe, not an actuator',
            retention='Captured flange, two M2x14 envelopes, side nuts. Sheath end crimp/bond, true threads and tools pending.'))
    parts['thorax_mount:reaction_frame']=frame
    candidate={**solids,sternum_key:sternum,**parts};changed={sternum_key,*parts};hits=[];inherited=[]
    for a,b in combinations(candidate,2):
        if a not in changed and b not in changed:continue
        volume=float(A.I.intersection_volume(candidate[a],candidate[b]))
        if volume<=.001:continue
        if sternum_key in (a,b):
            other=b if a==sternum_key else a
            baseline=float(A.I.intersection_volume(old,solids[other])) if other in solids else 0.
            new_overlap=float(A.I.intersection_volume(sternum-old,candidate[other]))
            if baseline>.001 and new_overlap<=.001:
                inherited.append(dict(other=other,source_overlap_mm3=baseline,current_overlap_mm3=volume,new_material_overlap_mm3=new_overlap));continue
        hits.append(dict(left=a,right=b,overlap_mm3=volume))
    probes=[]
    for row in ports:
        head=row['id'];wire=rod(row['wire_inboard_probe_source_mm'],row['ferrule_outlet_source_mm'],.3)
        for key,s in candidate.items():
            for kind,shape in [('wire_corridor',wire),('nominal_stub',stub_solids[head])]:
                v=float(A.I.intersection_volume(shape,s))
                if v>.001:probes.append(dict(id=head,kind=kind,target=key,overlap_mm3=v))
    components={k:sum(s.volume()>1e-8 for s in candidate[k].decompose()) for k in changed}
    report=dict(parent=provenance(parent_file),generator=provenance(ROOT/'tools/build_channel_first_thorax_reaction_mount.py'),
        mounts=mounts,ports=ports,source_fragment_repairs=source_repairs,new_interferences=hits,inherited_sternum_anatomy_junctions=inherited,
        probes=probes,positive_components=components,
        status='new_mount_static_screen_only_inherited_anatomy_not_mechanically_validated')
    if (dest/'fit_screen.json').exists():
        previous=(dest/'fit_screen.json').read_bytes()
        (dest/('fit_history_'+hashlib.sha256(previous).hexdigest()[:12]+'.json')).write_bytes(previous)
    save(dest/'fit_screen.json',report)
    if hits or probes or any(v!=1 for v in components.values()):
        print('Rejected chest mount; hits',hits,'probes',probes,'components',components,flush=True);return
    sources=dict(parent['sources']);origin=np.array([0.,-210.,1245.])
    for key in sorted(changed):
        _,source=export(candidate[key].as_original().simplify(1e-6),dest/(key.replace(':','_')+'.obj'),origin)
        sources[key]=source
    nominal_stubs={}
    for h,s in stub_solids.items():
        _,source=export(s,dest/(h+'_reference_stub.obj'),origin);nominal_stubs['thorax_stub:'+h]=source
    save(dest/'assembly.json',dict(sources=sources,ports=ports,mounts=mounts,
        nominal_conduit_stub_sources=nominal_stubs,arm_conduit_stub_sources=parent['conduit_stub_sources'],
        display_only_TPU_sources=parent['display_only_TPU_sources'],
        parent=provenance(parent_file),fit_screen=provenance(dest/'fit_screen.json'),
        generator=provenance(ROOT/'tools/build_channel_first_thorax_reaction_mount.py'),
        sternum_added_mm3=max(0.,float((sternum-old).volume())),sternum_removed_mm3=max(0.,float((old-sternum).volume())),
        status='local_sternum_two_port_candidate_not_motor_or_complete_supply',
        scope='No backboard. Original sternum kept except two nominal through-bores and local bosses; other bones/muscles unchanged. Inherited overlapping anatomical junction volumes explicitly retained, not whole-pair ignored. Full chest assembly, motor, sheath, tools, strength and base load path pending.'))
    print('Two sternum reaction ports exported;rigid parts',len(sources),'inherited junctions',len(inherited),flush=True)


if __name__=='__main__':main()
