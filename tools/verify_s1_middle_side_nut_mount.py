"""Independent middle nut slots, stock, contacts and conditional local assembly.

No complete muscle placement, physical threading or driven-ROM claim. A
clear local cap path may require subsequent wire feeding, recorded explicitly.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component,rod,shoulder_pads
import mesh_geom

DEST=A.DESIGN/'s1_middle_side_nut_mount'


def volume(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_middle_side_nut_mount.py')==r['script_sha256']
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
    assert sha(route_file)==r['routes_sha256']
    route=next(x for x in routes if x['head']=='FJ1467');port=route['ports']['origin']
    n=A.D.unit(port['normal_source']);I=np.array(port['inner_face_source_mm']);F=shoulder_pads.pad_basis(n)
    bone=solids['bone:scapula'];old_bone=load_component(parent['sources']['bone:scapula'])
    muscle=solids['tpu:deltoid'];cap_key='cap:FJ1467/origin';front_key='origin:FJ1467/front_support'
    stock=[]
    for i,hole in enumerate(r['selected_holes']):
        p=np.array(hole['axis_source_mm']);top=hole['nut_top_mm'];bottom=hole['nut_bottom_mm']
        roof=rod(p+n*(top+.3),p+n*(top+1.3),3.6)-rod(p+n*(top+.2),p+n*(top+1.4),1.8)
        floor=rod(p+n*(bottom-.9),p+n*(bottom-.3),3.6)-rod(p+n*(bottom-1.),p+n*(bottom-.2),1.8)
        nut=solids[f'origin:FJ1467/nut_{i}'];bolt=solids[f'origin:FJ1467/bolt_{i}']
        rotation=[]
        for angle in [-5.,5.]:rotation.append(dict(angle_deg=angle,overlap_mm3=volume(nut.transform(A.I.rotation_about(p,n,angle)),bone)))
        direction=np.array(r['nut_entries'][i]['selected']['direction_source'])
        stock.append(dict(index=i,roof_witness_missing_parent_mm3=float((roof-old_bone).volume()),
            roof_witness_missing_after_mm3=float((roof-bone).volume()),floor_witness_missing_after_mm3=float((floor-bone).volume()),
            nut_up_0p2mm_bone_overlap_mm3=volume(nut.translate(n*.2),bone),
            nut_down_0p2mm_bone_overlap_mm3=volume(nut.translate(-n*.2),bone),nut_rotation_probes=rotation,
            sideways_0p3mm_after_bolt_overlap_mm3=volume(nut.translate(direction*.3),bolt),
            nut_entry_start_bone_overlap_mm3=volume(nut.translate(direction*r['nut_entries'][i]['selected']['corridor_mm']),bone)))
    stock_pass=all(x['roof_witness_missing_after_mm3']<.001 and x['floor_witness_missing_after_mm3']<.001
        and x['nut_up_0p2mm_bone_overlap_mm3']>.001 and x['nut_down_0p2mm_bone_overlap_mm3']>.001
        and all(y['overlap_mm3']>.001 for y in x['nut_rotation_probes'])
        and x['sideways_0p3mm_after_bolt_overlap_mm3']>.001 and x['nut_entry_start_bone_overlap_mm3']<.001 for x in stock)
    contacts=dict(front_bone_mm3=volume(solids[front_key].translate(-n*.02),bone),
        front_TPU_mm3=volume(solids[front_key].translate(n*.02),muscle),
        cap_TPU_mm3=volume(solids[cap_key].translate(-n*.1),muscle),
        cap_TPU_gap_capped_2mm=float(A.I.gap(solids[cap_key],muscle,2.)))
    contact_pass=all(contacts[k]>.001 for k in ['front_bone_mm3','front_TPU_mm3','cap_TPU_mm3'])
    fixed=set(solids)-set(r['new_keys'])-{'tpu:deltoid',cap_key};stages=[]
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
    for i,entry in enumerate(r['nut_entries']):
        replay(f'side_nut_{i}_before_muscle',[f'origin:FJ1467/nut_{i}'],
            np.array(entry['selected']['direction_source']),entry['selected']['corridor_mm'],.1)
    replay('front_support_before_muscle',[front_key],n,15.)
    fixed.add('tpu:deltoid')  # Actual muscle insertion is explicitly unverified.
    search=[];cap=solids[cap_key]
    for elevation in [75.,60.,45.,30.,15.,0.]:
        for azimuth in range(0,360,15):
            a=np.deg2rad(azimuth);b=np.deg2rad(elevation)
            direction=(F[:,0]*np.cos(a)+F[:,1]*np.sin(a))*np.cos(b)+n*np.sin(b)
            hit=None
            for s in np.arange(0.,12.0001,.5):
                posed=cap.translate(direction*s)
                for key in sorted(fixed):
                    v=volume(posed,solids[key])
                    if v>.001:hit=dict(offset_mm=float(s),target=key,volume_mm3=v);break
                if hit:break
            search.append(dict(elevation_deg=elevation,azimuth_deg=azimuth,direction_source=direction.tolist(),first_collision=hit,clear=hit is None))
        print('cap elevation',elevation,'clear',sum(x['clear'] for x in search if x['elevation_deg']==elevation),flush=True)
    eligible=[x for x in search if x['clear']];selected=eligible[0] if eligible else None
    if selected:replay('middle_cap_before_wire',[cap_key],np.array(selected['direction_source']),12.,.1)
    else:fixed.add(cap_key)  # Explicit failed placement, later checks conditional.
    for i in range(2):
        replay(f'front_screw_{i}',[f'origin:FJ1467/bolt_{i}',f'origin:FJ1467/front_washer_{i}'],n,30.)
    assert fixed==set(solids)
    wire_hits=[]
    if selected:
        direction=np.array(selected['direction_source'])
        for route in routes:
            corridor=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3);hits=[]
            for s in np.arange(0.,12.0001,.1):
                v=volume(cap.translate(direction*s),corridor)
                if v>.001:hits.append(dict(offset_mm=float(s),volume_mm3=v))
            if hits:wire_hits.append(dict(head=route['head'],samples=hits))
    tools=[];thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    for i,h in enumerate(r['selected_holes']):
        p=np.array(h['axis_source_mm']);tool=rod(p+n*(thickness+3.9),p+n*(thickness+25.),1.5)
        for key,solid in solids.items():
            v=volume(tool,solid)
            if v>.001:tools.append(dict(index=i,target=key,volume_mm3=v))
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),topology=topology,
        topology_pass=topology_pass,unchanged_parent_rows=len(unchanged),post_slot_stock_and_capture=stock,
        nominal_nut_roof_floor_antirotation_pass=stock_pass,contact_probe_shifts_mm=dict(front=.02,cap=.1),contact_probes=contacts,
        nominal_contact_chain_present=contact_pass,cap_direction_search=search,selected_cap_entry=selected,
        cap_direction_step_deg=15.,cap_coarse_step_mm=.5,cap_must_precede_wires=[x['head'] for x in wire_hits],
        cap_entry_wire_corridor_hits=wire_hits,local_stages=stages,poses=sum(x['poses'] for x in stages),all_objects_accounted_for=True,
        local_approaches_geometry_pass=bool(selected) and not any(x['hits'] for x in stages),
        tool_shaft_envelope_hits=tools,real_tool_engagement_verified=False,
        static_fit_stock_contact_pass=r['static_geometry_pass'] and topology_pass and stock_pass and contact_pass,
        whole_module_installation_verified=False,muscle_placement_transition_verified=False,physical_wire_feeding_verified=False,
        manual_nut_retention_before_bolts_verified=False,real_threads_preload_strength_verified=False,
        driven_motion_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print('stock',stock,'contacts',contacts,'static',result['static_fit_stock_contact_pass'],'paths',result['local_approaches_geometry_pass'],flush=True)


if __name__=='__main__':main()
