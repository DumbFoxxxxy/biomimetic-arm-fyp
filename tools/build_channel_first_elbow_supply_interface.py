"""Positive ferrule capture on original elbow caps, without relocating muscles."""
from itertools import combinations
import argparse
import hashlib
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,export,rod
from build_channel_first_inventory import provenance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flexor-ear-deg',type=float,default=0.)
    parser.add_argument('--extensor-ear-deg',type=float,default=0.)
    parser.add_argument('--ear-layout',choices=['symmetric','free_side'],default='free_side')
    args=parser.parse_args()
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'elbow_extension_stop/assembly.json';parent=read(parent_file)
    mount_file=base/'elbow_mount_candidate/assembly.json';mount=read(mount_file)
    dest=base/'elbow_supply_interface';dest.mkdir(exist_ok=True)
    original={k:load_component(v) for k,v in parent['sources'].items()}
    candidate=dict(original);changed=set();rows=[];sheaths={};mf=A.I._manifold()
    for muscle in mount['muscles'].values():
        head=muscle['anatomical_id'];path=np.array(muscle['wire_path_source_mm'])
        P=path[1];axis=path[0]-P;axis/=np.linalg.norm(axis)
        design=read(ROOT/muscle['source_definition']['file'])
        x=np.array(design['basis_columns_display'])[:,0];assert abs(x@axis)<1e-9
        y=np.cross(axis,x)
        angle=np.radians(args.flexor_ear_deg if head=='FJ1486' else args.extensor_ear_deg)
        x,y=x*np.cos(angle)+y*np.sin(angle),-x*np.sin(angle)+y*np.cos(angle)
        F=np.stack([x,y,axis],axis=1)
        def pt(a,b,t):return P+F@np.array([a,b,t])
        def cyl(a,b,t0,t1,r):return rod(pt(a,b,t0),pt(a,b,t1),r)
        key=f'elbow:{head}_proximal_rigid_cap';old=original[key]
        stock=cyl(0,0,20.65,32.35,4.)
        cuts=[cyl(0,0,20.55,25.15,.85),cyl(0,0,25.15,31.15,2.1),cyl(0,0,31.15,32.45,3.1)]
        plate=cyl(0,0,32.35,34.35,4.)-cyl(0,0,32.25,34.45,1.7)
        additions={}
        if args.ear_layout=='free_side':
            free_sign=-1. if head=='FJ1486' else 1.
            positions=[(s*9/np.sqrt(2),free_sign*9/np.sqrt(2)) for s in (-1.,1.)]
        else:positions=[(0.,s*6.) for s in (-1.,1.)]
        for i,(a,b) in enumerate(positions):
            radius=np.hypot(a,b);normal=np.array([a,b])/radius;tangent=np.array([normal[1],-normal[0]])
            H=F@np.array([[tangent[0],normal[0],0],[tangent[1],normal[1],0],[0,0,1]])
            def earbox(lo,hi):return mf.Manifold.cube(np.array(hi)-lo).translate(lo).transform(np.c_[H,pt(a,b,0)])
            def earhex(t0,height,af):return mf.Manifold.cylinder(height,af/np.sqrt(3),af/np.sqrt(3),6).rotate([0,0,30]).transform(np.c_[H,pt(a,b,t0)])
            stock+=cyl(a,b,20.65,32.35,3.2)
            if radius>7.:
                stock+=rod(pt(*(normal*3.),27.),pt(a,b,27.),1.6)
                plate+=rod(pt(*(normal*3.),33.35),pt(a,b,33.35),1.)
            # Side-loadable nut pocket; the M2 bolt then prevents radial escape.
            # One rectangular side slot avoids coincident stepped pocket cuts.
            # Its4.2mm tangential width captures the4mm nut flats; nut bears
            # directly on the22.9mm ceiling, without a buried loose washer.
            pocket=earbox([-2.1,-2.6,20.9],[2.1,5.,22.9])
            cuts.append(pocket)
            hole=cyl(a,b,20.55,37.,1.15)
            cuts.append(hole)
            plate+=(cyl(a,b,32.35,34.35,3.2)-hole)
            plate-=hole
            stem=f'supply:{head}_'
            additions[stem+f'bolt_{i}']=cyl(a,b,20.65,34.65,1.)+cyl(a,b,34.65,36.65,1.9)
            additions[stem+f'head_washer_{i}']=cyl(a,b,34.35,34.65,2.5)-hole
            additions[stem+f'nut_{i}']=earhex(21.3,1.6,4.)-cyl(a,b,21.2,23.,1.05)
        # Retain the old0.8mm bore through original cap; the new0.85mm throat
        # cuts the added boss only. Two M2 bores reach0.5mm into its top face.
        centre_cuts=mf.Manifold.batch_boolean(cuts[:3],mf.OpType.Add)
        local_cuts=mf.Manifold.batch_boolean(cuts[3:],mf.OpType.Add)
        new_cap=(old+(stock-centre_cuts))-local_cuts
        ferrule=(cyl(0,0,25.15,31.15,2.)+cyl(0,0,31.15,32.15,3.))
        ferrule-=cyl(0,0,25.05,32.25,.65)+cyl(0,0,28.15,32.25,1.55)
        additions[f'supply:{head}_ferrule']=ferrule
        additions[f'supply:{head}_retainer']=plate
        candidate[key]=new_cap;candidate.update(additions);changed.update([key,*additions])
        sheaths[head]=cyl(0,0,28.15,45.15,1.5)-cyl(0,0,28.05,45.25,.6)
        rows.append(dict(id=head,cap_key=key,frame_columns=F.tolist(),origin_source_mm=P.tolist(),
            original_wire_tail_source_mm=path[0].tolist(),reaction_axis_source=axis.tolist(),
            ferrule_bearing_point_source_mm=pt(0,0,25.15).tolist(),
            sheath_bearing_point_source_mm=pt(0,0,28.15).tolist(),
            supply_exit_source_mm=pt(0,0,45.15).tolist(),supply_exit_tangent_outward=axis.tolist(),
            old_cap_removed_mm3=max(0.,float((old-new_cap).volume())),added_cap_mm3=max(0.,float((new_cap-old).volume())),
            ear_centres_local_mm=[list(c) for c in positions],
            ferrule_seat_radius_mm=2.1,ferrule_outer_radius_mm=2.,flange_radius_mm=3.,
            maximum_ferrule_axial_float_mm=.2,sheath_outer_radius_mm=1.5,sheath_lumen_radius_mm=.6,
            wire_radius_mm=.3,bearing_annulus_mm2=float(np.pi*(2.**2-.85**2)),
            owner='humerus',
            reaction_chain='sheath end -> ferrule step -> integral proximal-cap boss -> original cap M3s -> original bone frame -> humerus M3 mount seats',
            fixation='ferrule flange captured by removable plate and two M2x14 envelope bolts, captive side-loaded M2 nuts',
            unresolved='Sheath-to-ferrule crimp/bond retention, real thread engagement, strength, tolerances and assembly access not yet verified. Sheath stub material/compression law unspecified; not a bare TPU compression-proof claim.'))
    hits=[]
    for a,b in combinations(candidate,2):
        if a not in changed and b not in changed:continue
        v=float(A.I.intersection_volume(candidate[a],candidate[b]))
        if v>.001:hits.append(dict(left=a,right=b,overlap_mm3=v))
    components={k:sum(c.volume()>1e-8 for c in candidate[k].decompose()) for k in changed}
    flexible_hits=[];wire_hits=[]
    for r in rows:
        h=r['id'];P=np.array(r['origin_source_mm']);axis=np.array(r['reaction_axis_source'])
        wire=rod(P-axis*.1,P+axis*45.15,.3)
        for k,s in candidate.items():
            v=float(A.I.intersection_volume(sheaths[h],s))
            if v>.001:flexible_hits.append(dict(channel=h,target=k,overlap_mm3=v))
            v=float(A.I.intersection_volume(wire,s))
            if v>.001:wire_hits.append(dict(channel=h,target=k,overlap_mm3=v))
    report=dict(parent=provenance(parent_file),muscle_reference=provenance(mount_file),parameters=vars(args),
        generator=provenance(ROOT/'tools/build_channel_first_elbow_supply_interface.py'),
        interfaces=rows,changed_components=components,rigid_interferences=hits,
        sheath_stub_interferences=flexible_hits,wire_entry_interferences=wire_hits)
    if (dest/'fit_screen.json').exists():
        previous=(dest/'fit_screen.json').read_bytes()
        (dest/('fit_history_'+hashlib.sha256(previous).hexdigest()[:12]+'.json')).write_bytes(previous)
    save(dest/'fit_screen.json',report)
    if hits or wire_hits or flexible_hits or any(v!=1 for v in components.values()):
        print('Rejected inlet candidate; rigid',hits,'stub',flexible_hits,'steel',wire_hits,'components',components,flush=True);return
    sources=dict(parent['sources'])
    for key in sorted(changed):
        _,meta=export(candidate[key].as_original().simplify(1e-6),dest/(key.replace(':','_')+'.obj'),np.mean([r['origin_source_mm'] for r in rows],axis=0))
        sources[key]=meta
    conduit_sources={}
    for h,s in sheaths.items():
        _,meta=export(s,dest/(h+'_sheath_stub.obj'),np.array(next(r['origin_source_mm'] for r in rows if r['id']==h)))
        conduit_sources['sheath:'+h]=meta
    reloaded={k:load_component(s) for k,s in sources.items()};roundtrip_hits=[]
    for a,b in combinations(reloaded,2):
        if a not in changed and b not in changed:continue
        v=float(A.I.intersection_volume(reloaded[a],reloaded[b]))
        if v>.001:roundtrip_hits.append(dict(left=a,right=b,overlap_mm3=v))
    for r in rows:
        P=np.array(r['origin_source_mm']);axis=np.array(r['reaction_axis_source'])
        wire=rod(P-axis*.1,P+axis*45.15,.3);stub=load_component(conduit_sources['sheath:'+r['id']])
        for k,s in reloaded.items():
            for label,shape in [('steel',wire),('sheath_stub',stub)]:
                v=float(A.I.intersection_volume(shape,s))
                if v>.001:roundtrip_hits.append(dict(channel=r['id'],kind=label,target=k,overlap_mm3=v))
    save(dest/'export_screen.json',dict(sources_checked=len(sources),interferences=roundtrip_hits,
        source_fingerprints={k:v['sha256'] for k,v in sources.items()},scope='Reimported exact208 rigid meshes and two nominal stub tubes; same neutral threshold0.001mm3. No source-pair exclusions.'))
    assert not roundtrip_hits,roundtrip_hits
    save(dest/'assembly.json',dict(sources=sources,conduit_stub_sources=conduit_sources,
        display_only_TPU_sources=parent['display_only_TPU_sources'],interfaces=rows,
        parent=provenance(parent_file),fit_screen=provenance(dest/'fit_screen.json'),
        export_screen=provenance(dest/'export_screen.json'),
        generator=provenance(ROOT/'tools/build_channel_first_elbow_supply_interface.py'),
        status='neutral_supply_reaction_interface_candidate_full_route_assembly_and_motion_pending',
        scope='Rigid ferrule capture and end reaction only. No chest motor, full supply sheath or claimed flexible dynamics. Original muscles, local wire path and bone geometry retained. End crimp/retention and nominal fastener threads need physical design/validation.'))
    save(dest/'view_scene.json',dict(sources={**sources,**conduit_sources,**parent['display_only_TPU_sources']},
        status='static_supply_inlet_candidate_not_full_route'))
    print('Two reaction/capture interfaces exported, rigid parts',len(sources),flush=True)


if __name__=='__main__':main()
