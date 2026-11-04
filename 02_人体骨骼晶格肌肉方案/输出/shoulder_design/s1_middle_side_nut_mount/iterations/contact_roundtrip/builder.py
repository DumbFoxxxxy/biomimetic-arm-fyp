"""Middle deltoid clamp with actual side-loaded nut slots in the scapula.

Keeps bone/muscle placements and the original wire port. Relieves only the
declared terminal to fit a real support and flange. Slots are cut, never
excused from collision tests. Strength, threads and driven motion unverified.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,export,load_component,shoulder_pads

DEST=A.DESIGN/'s1_middle_side_nut_mount'


def main():
    DEST.mkdir(exist_ok=True)
    study_file=A.DESIGN/'s1_middle_blind_study/side_nut_study.json';study=read(study_file)
    assert study['script_sha256']==sha(ROOT/'tools/study_s1_middle_side_nuts.py')
    choice=study['selected'];assert choice and study['screen_pass']
    file=ROOT/study['parent_report'];parent=read(file);assert sha(file)==study['parent_sha256']
    original={k:load_component(s) for k,s in parent['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==study['routes_sha256']
    port=next(r for r in routes if r['head']=='FJ1467')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n);mf=A.I._manifold()
    bone=original['bone:scapula'];muscle=original['tpu:deltoid'];old_cap=original['cap:FJ1467/origin']
    zone=rod(I-n*18.,I+n*.2,12.);cap=old_cap
    for sign in [-1,1]:
        p=I+F[:,0]*sign*7.
        cap=cap+(rod(p-n*.05,p+n*(thickness+.05),1.71)^rod(I,I+n*thickness,11.))
    parts={};drilled=bone;entries=[];pockets=[]
    # Cut side loading before muscles/supports are fitted. Final slots are
    # retained in the source scapula, and checked again after export.
    early_obstacles={k:s for k,s in original.items() if k not in ['bone:scapula','tpu:deltoid','cap:FJ1467/origin']}
    for index,hole in enumerate(choice['holes']):
        p=np.array(hole['axis_source_mm']);bottom=hole['nut_bottom_mm'];top=hole['nut_top_mm']
        hardware=M.m3_parts(p,n,thickness,bottom+3.15,hole['screw_length_mm'])
        nut=hardware['nut'];candidates=[]
        for angle in range(0,360,45):
            a=np.deg2rad(angle);direction=F[:,0]*np.cos(a)+F[:,1]*np.sin(a)
            distance=next((float(s) for s in np.arange(8.,25.001,.5)
                           if A.I.intersection_volume(nut.translate(direction*s),bone)<.001),None)
            if distance is None:continue
            hits=[]
            for s in np.arange(0.,distance+.001,.5):
                posed=nut.translate(direction*s)
                for key,obstacle in early_obstacles.items():
                    v=float(A.I.intersection_volume(posed,obstacle))
                    if v>.001:hits.append(dict(offset_mm=float(s),target=key,volume_mm3=v))
            candidates.append(dict(azimuth_deg=angle,direction_source=direction.tolist(),corridor_mm=distance,hits=hits))
        clear=[x for x in candidates if not x['hits']]
        if not clear:
            save(DEST/'entry_failure.json',dict(index=index,candidates=candidates,stage_pass=False))
            raise ValueError('No side nut entry before deltoid; saved evidence')
        selected=min(clear,key=lambda x:x['corridor_mm']);direction=np.array(selected['direction_source'])
        pocket=mf.Manifold.cylinder(2.7,5.7/np.sqrt(3),5.7/np.sqrt(3),6).translate([0,0,bottom-.15]).transform(np.c_[F,p])
        sweep=mf.Manifold.batch_hull([pocket,pocket.translate(direction*(selected['corridor_mm']+.5))])
        drilled=drilled-sweep;pockets.append(sweep)
        # Blind bolt continuation gives a nominal .5mm tip reserve while
        # retaining the measured deep bone face and GH clearance.
        bore=rod(p+n*(hole['screw_tip_mm']-.5),p+n*(thickness+5.),1.7)
        drilled=drilled-bore
        for name in ['nut','bolt','front_washer']:parts[f'origin:FJ1467/{name}_{index}']=hardware[name]
        cap=cap+rod(p,p+n*thickness,4.)
        cap=cap-bore
        entries.append(dict(index=index,axis_source_mm=p.tolist(),nut_bottom_mm=bottom,nut_top_mm=top,
            selected=selected,candidates=candidates,screw_length_mm=hole['screw_length_mm'],screw_tip_mm=hole['screw_tip_mm']))
        print('nut',index,'entry',selected['azimuth_deg'],selected['corridor_mm'],flush=True)
    # Bone is shallower here than at anterior/posterior origins. A copied
    # deep spacer would enter the scapula. Use a shallow conforming support.
    front=rod(I-n*3.5,I-n*1.5,10.8)
    for hole in choice['holes']:
        p=np.array(hole['axis_source_mm']);face=hole['ray_crossings_normal_mm'][0]
        front=front+rod(p+n*(face-.2),p-n*1.5,3.8)
    front=front-bone-rod(I-n*18.,I+n*.2,3.)
    flange=rod(I-n*1.5,I-n*.03,10.3)
    for r in routes:
        path=np.array(r['full_wire_path_source_mm'])
        front=front-A.D.outer_envelope(path,2.1)
        flange=flange-A.D.outer_envelope(path,.86)
        cap=cap-A.D.outer_envelope(path,.86)
    revised=(muscle+flange)-(front^zone)-(bone^zone)
    cap=cap-revised
    for hole in choice['holes']:
        p=np.array(hole['axis_source_mm']);drill=rod(p-n*20.,p+n*10.,1.7)
        front=front-drill;revised=revised-(drill^zone)
    parts.update({'origin:FJ1467/front_support':front,'bone:scapula':drilled,
                  'tpu:deltoid':revised,'cap:FJ1467/origin':cap})
    sources=dict(parent['sources']);actual=dict(original);cleanup=[]
    for key,solid in parts.items():
        prepared=solid.as_original().simplify(1e-6)
        cleanup.append(dict(key=key,added_mm3=float((prepared-solid).volume()),removed_mm3=float((solid-prepared).volume())))
        origin=np.array(sources[key]['translation_source_mm'] or [0.,0.,0.]) if key in sources else I
        actual[key],sources[key]=export(prepared,DEST/(key.replace(':','_').replace('/','_')+'.obj'),origin)
        print('exported',key,flush=True)
    pairs=[];count=0;keys=list(actual)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in parts and b not in parts:continue
            count+=1;v=float(A.I.intersection_volume(actual[a],actual[b]))
            previous=float(A.I.intersection_volume(original[a],original[b])) if a in original and b in original else 0.
            if v>.001 or previous>.001:pairs.append(dict(a=a,b=b,volume_mm3=v,parent_mm3=previous,new_or_increased=v>max(.001,previous+.001)))
    wires=[]
    for r in routes:
        corridor=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3);hits=[]
        for key in parts:
            v=float(A.I.intersection_volume(actual[key],corridor))
            if v>.001:hits.append(dict(key=key,volume_mm3=v))
        wires.append(dict(head=r['head'],hits=hits))
    removed=muscle-actual['tpu:deltoid'];own=next(r for r in routes if r['head']=='FJ1467')
    wall_loss=float(A.I.intersection_volume(removed,A.D.outer_envelope(np.array(own['full_wire_path_source_mm']),2.)))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),study_sha256=sha(study_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(set(parts)-set(original)),replaced_keys=sorted(set(parts)&set(original)),
        unchanged_parent_rows=len(set(original)-set(parts)),cap_centre_and_wire_unchanged=True,
        selected_holes=choice['holes'],nut_entries=entries,numerical_cleanup=cleanup,
        dimensions=dict(pitch_mm=choice['pitch_mm'],nut_AF_mm=5.5,nut_height_mm=2.4,pocket_AF_mm=5.7,
            pocket_height_mm=2.7,nominal_roof_stock_mm=1.5,nominal_floor_stock_mm=1.,
            front_support_normal_mm=[-3.5,-1.5],flange_normal_mm=[-1.5,-.03],
            terminal_zone_radius_mm=12.,terminal_zone_normal_mm=[-18.,.2]),
        bone_removed_mm3=float((bone-actual['bone:scapula']).volume()),
        TPU_added_mm3=float((actual['tpu:deltoid']-muscle).volume()),TPU_removed_mm3=float(removed.volume()),
        TPU_removed_outside_terminal_mm3=float((removed-zone).volume()),working_wire_R2_wall_removed_mm3=wall_loss,
        affected_pairs=count,material_intersections=pairs,six_active_wire_corridors=wires,
        static_geometry_pass=not any(x['new_or_increased'] for x in pairs) and not any(x['hits'] for x in wires) and wall_loss<.001,
        intended_body_owners={**parent['intended_body_owners'],**{k:'scapula_ideal_fastened' for k in parts if k.startswith('origin:')}},
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None,
        not_verified=['post-cut roof and floor witnesses','full staged assembly','real threads, preload and strength',
                      'manual nut retention until bolt fits','whole muscle placement, wire feeding and actuation'])
    save(DEST/'assembly_check.json',report)
    print('instances',len(sources),'conflicts',[x for x in pairs if x['new_or_increased']],'wire wall loss',wall_loss,flush=True)


if __name__=='__main__':main()
