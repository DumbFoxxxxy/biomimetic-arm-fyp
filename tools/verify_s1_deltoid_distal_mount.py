"""Shared distal clamp: readback, bearing, tools and conditional assembly.

Three historical stop blanks are explicitly unverified references, placed
after the cap. They are retained in the final scene, never counted as real
wire stops or validated installation. Muscle insertion is also unverified.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component,rod,shoulder_pads
import mesh_geom

DEST=A.DESIGN/'s1_deltoid_distal_mount'


def volume(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_deltoid_distal_mount.py')==r['script_sha256']
    changed=set(r['new_keys'])|set(r['replaced_keys']);unchanged=set(parent['sources'])-changed
    for key in unchanged:assert r['sources'][key]==parent['sources'][key]
    solids={k:load_component(s) for k,s in r['sources'].items()};topology=[]
    for key in sorted(changed):
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']))
        _,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f];_,boundary,nm=mesh_geom.edge_census(faces)
        xyz=v[f];repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        zero=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        count=sum(p.volume()>1e-8 for p in solids[key].decompose())
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nm,repeated=repeated,zero_area=zero,positive_components=count))
    topology_pass=all(x['boundary']==x['nonmanifold']==x['repeated']==x['zero_area']==0 and x['positive_components']==1 for x in topology)
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==r['routes_sha256'];port=next(x for x in routes if x['head']=='FJ1468')['ports']['distal']
    n=A.D.unit(port['normal_source']);I=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm']);F=shoulder_pads.pad_basis(n)
    mf=A.I._manifold();front_key='distal:FJ1468/front_support';rear_key='distal:FJ1468/backing_plate';cap_key='cap:FJ1468/distal'
    bone=solids['bone:humerus'];muscle=solids['tpu:deltoid'];front=solids[front_key];back=solids[rear_key];cap=solids[cap_key]
    contacts=dict(front_bone_mm3=volume(front.translate(-n*.02),bone),front_TPU_mm3=volume(front.translate(n*.02),muscle),
        rear_bone_mm3=volume(back.translate(n*.02),bone),cap_TPU_mm3=volume(cap.translate(-n*.1),muscle),
        cap_TPU_gap_capped_2mm=float(A.I.gap(cap,muscle,2.)))
    contact_pass=all(contacts[k]>.001 for k in ['front_bone_mm3','front_TPU_mm3','rear_bone_mm3','cap_TPU_mm3'])
    bearing=[]
    for i,hole in enumerate(r['selected_holes']):
        bolt=solids[f'distal:FJ1468/bolt_{i}'];fw=solids[f'distal:FJ1468/front_washer_{i}'];bw=solids[f'distal:FJ1468/back_washer_{i}'];nut=solids[f'distal:FJ1468/nut_{i}']
        bearing.append(dict(index=i,front_head_to_washer_mm3=volume(bolt.translate(-n*.12),fw),
            front_washer_to_cap_mm3=volume(fw.translate(-n*.2),cap),rear_nut_to_washer_mm3=volume(nut.translate(n*.2),bw),
            rear_washer_to_backing_mm3=volume(bw.translate(n*.2),back)))
    bearing_pass=all(all(v>.001 for k,v in row.items() if k!='index') for row in bearing)
    terminal_holding=[]
    for direction in [F[:,0],-F[:,0],F[:,1],-F[:,1]]:
        posed=muscle.translate(direction*.3)
        terminal_holding.append(dict(direction_source=direction.tolist(),
            overlaps=[dict(key=key,volume_mm3=volume(posed,solids[key])) for key in [front_key,cap_key,'distal:FJ1468/bolt_0','distal:FJ1468/bolt_1']]))
    holding_pass=all(any(h['volume_mm3']>.001 for h in x['overlaps']) for x in terminal_holding)
    blanks={f'stop:{head}_stop_blank' for head in ['FJ1468','FJ1513','FJ1467']};assert blanks<=set(solids)
    stop_seats=[]
    for key in sorted(blanks):
        sleeve=solids[key]
        stop_seats.append(dict(key=key,neutral_cap_overlap_mm3=volume(sleeve,cap),
            inward_0p06mm_cap_overlap_mm3=volume(sleeve.translate(-n*.06),cap)))
    stop_seat_pass=all(x['neutral_cap_overlap_mm3']<.001 and x['inward_0p06mm_cap_overlap_mm3']>.001 for x in stop_seats)
    fixed=set(solids)-set(r['new_keys'])-{'tpu:deltoid',cap_key}-blanks;stages=[]
    def replay(name,keys,direction,start,step=.5):
        assert not set(keys)&fixed;hits=[];offsets=np.arange(0.,start+.0001,step)
        for way,values in [('in',offsets[::-1]),('out',offsets)]:
            for s in values:
                for key in keys:
                    posed=solids[key].translate(direction*s)
                    for target in sorted(fixed):
                        v=volume(posed,solids[target])
                        if v>.001:hits.append(dict(way=way,offset_mm=float(s),moving=key,target=target,volume_mm3=v))
        stages.append(dict(name=name,moving=keys,already_installed=sorted(fixed),direction_source=direction.tolist(),
            step_mm=step,poses=2*len(offsets),hits=hits,conditional_on_prior_stages=True))
        fixed.update(keys);print(name,'hits',len(hits),flush=True)
    replay('front_support_before_deltoid',[front_key],n,20.)
    fixed.add('tpu:deltoid')  # Actual muscle insertion is not a checked stage.
    search=[]
    for elevation in [75.,60.,45.,30.,15.,0.]:
        for azimuth in range(0,360,15):
            a=np.deg2rad(azimuth);b=np.deg2rad(elevation)
            direction=(F[:,0]*np.cos(a)+F[:,1]*np.sin(a))*np.cos(b)+n*np.sin(b);hit=None
            for s in np.arange(0.,16.0001,.5):
                posed=cap.translate(direction*s)
                for key in sorted(fixed):
                    v=volume(posed,solids[key])
                    if v>.001:hit=dict(offset_mm=float(s),target=key,volume_mm3=v);break
                if hit:break
            search.append(dict(elevation_deg=elevation,azimuth_deg=azimuth,direction_source=direction.tolist(),first_collision=hit,clear=hit is None))
        print('cap elevation',elevation,'clear',sum(x['clear'] for x in search if x['elevation_deg']==elevation),flush=True)
    eligible=[x for x in search if x['clear']];selected=eligible[0] if eligible else None
    if selected:replay('shared_cap_before_wires_and_stops',[cap_key],np.array(selected['direction_source']),16.,.1)
    else:fixed.add(cap_key)
    replay('rear_backing',[rear_key],-n,35.)
    for i in range(2):replay(f'front_screw_{i}',[f'distal:FJ1468/bolt_{i}',f'distal:FJ1468/front_washer_{i}'],n,45.)
    for i in range(2):replay(f'rear_nut_{i}',[f'distal:FJ1468/back_washer_{i}',f'distal:FJ1468/nut_{i}'],-n,15.)
    fixed.update(blanks)  # Unverified stop blank transition; not a placement pass.
    assert fixed==set(solids)
    wire_hits=[]
    if selected:
        direction=np.array(selected['direction_source'])
        for route in routes:
            corridor=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3);hits=[]
            for s in np.arange(0.,16.0001,.1):
                v=volume(cap.translate(direction*s),corridor)
                if v>.001:hits.append(dict(offset_mm=float(s),volume_mm3=v))
            if hits:wire_hits.append(dict(head=route['head'],samples=hits))
    tool_hits=[];th=r['dimensions']['cap_thickness_mm'];lower=r['dimensions']['rear_plate_lower_mm']
    for i,hole in enumerate(r['selected_holes']):
        p=np.array(hole['axis_source_mm']);front_tool=rod(p+n*(th+3.9),p+n*(th+25.),1.5)
        top=lower-.75-.1;base=top-8.
        socket=rod(p+n*base,p+n*top,4.)
        cavity=mf.Manifold.cylinder(8.2,5.7/np.sqrt(3),5.7/np.sqrt(3),6).translate([0,0,base-.1]).transform(np.c_[F,p])
        rear_tool=(socket-cavity)+rod(p+n*(base-25.),p+n*(base+.5),2.)
        for name,tool in [('front_driver_shaft',front_tool),('rear_nominal_socket_OD8',rear_tool)]:
            for key,solid in solids.items():
                v=volume(tool,solid)
                if v>.001:tool_hits.append(dict(index=i,tool=name,target=key,volume_mm3=v))
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),topology=topology,topology_pass=topology_pass,
        unchanged_parent_rows=len(unchanged),contact_probe_shifts_mm=dict(front=.02,rear=.02,cap=.1),contact_probes=contacts,
        nominal_contact_chain_present=contact_pass,fastener_bearing_probes=bearing,fastener_bearing_pass=bearing_pass,
        terminal_shear_probes=terminal_holding,nominal_terminal_shear_restraint_present=holding_pass,
        original_sleeve_seat_probes=stop_seats,original_sleeve_seat_geometry_preserved=stop_seat_pass,
        cap_direction_search=search,selected_cap_entry=selected,cap_direction_step_deg=15.,cap_coarse_step_mm=.5,
        cap_must_precede_wires=[x['head'] for x in wire_hits],cap_entry_wire_corridor_hits=wire_hits,
        local_stages=stages,poses=sum(x['poses'] for x in stages),all_objects_accounted_for=True,
        unverified_late_stop_blank_keys=sorted(blanks),stop_blanks_counted_as_real_terminations=False,
        local_approaches_geometry_pass=bool(selected) and not any(x['hits'] for x in stages),
        tool_envelope_hits=tool_hits,tool_envelope_pass=not tool_hits,real_tool_engagement_verified=False,
        static_fit_contact_bearing_pass=r['static_geometry_pass'] and topology_pass and contact_pass and bearing_pass and holding_pass and stop_seat_pass,
        whole_module_installation_verified=False,muscle_placement_transition_verified=False,physical_wire_feeding_verified=False,
        real_threads_preload_strength_verified=False,driven_motion_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('contacts',contacts,'bearing',bearing,'static',result['static_fit_contact_bearing_pass'],
          'paths',result['local_approaches_geometry_pass'],'tools',tool_hits,flush=True)


if __name__=='__main__':main()
