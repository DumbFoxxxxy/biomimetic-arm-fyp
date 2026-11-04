"""Find two fully supported blind inserts inside the middle terminal area.

The symmetric pitch-14 sweep had no adequate radial bone stock. This study
keeps the muscle/cap centre and port, but allows the two fastener positions
to move independently in a 2mm grid. No source geometry is changed.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_middle_blind_study'


def main():
    first_file=DEST/'blind_insert_study.json';first=read(first_file)
    assert not first['screen_pass']
    file=ROOT/first['parent_report'];parent=read(file);assert sha(file)==first['parent_sha256']
    solids={k:load_component(s) for k,s in parent['sources'].items()};bone=solids['bone:scapula']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==first['routes_sha256']
    port=next(r for r in routes if r['head']=='FJ1467')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    zone=rod(I-n*18.,I+n*.2,12.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    rows=[]
    for x in range(-10,11,2):
        for y in range(-10,11,2):
            if x*x+y*y>100:continue
            p=I+F[:,0]*x+F[:,1]*y;crossings=bone.ray_cast(p+n*10.,p-n*80.)
            zs=[float((np.array(h.position)-I)@n) for h in crossings]
            for setback in [1.,2.,3.,4.]:
                hits=[]
                row=dict(offset_local_mm=[x,y],axis_source_mm=p.tolist(),setback_mm=setback,ray_crossings_normal_mm=zs)
                if len(zs)<2:
                    rows.append(dict(**row,hits=[dict(reason='no first segment')],preliminary_clear=False));continue
                top=zs[0]-setback;bottom=top-6.;reserve=bottom-.5-zs[1]
                guard=rod(p+n*bottom,p+n*top,3.5);missing=float((guard-bone).volume())
                row.update(insert_top_mm=top,insert_bottom_mm=bottom,radial_guard_missing_mm3=missing,bottom_axial_reserve_mm=reserve)
                if missing>.001:hits.append(dict(reason='radial/bottom guard not contained',missing_mm3=missing))
                if reserve<1.:hits.append(dict(reason='less than 1mm axial reserve',reserve_mm=reserve))
                if hits:
                    rows.append(dict(**row,hits=hits,preliminary_clear=False));continue
                options=[]
                for length in [10.,12.,14.,16.,18.,20.,22.,25.]:
                    tip=thickness+.75-length;engaged=top-max(bottom,tip)
                    if engaged>=4. and tip>=bottom-.3:options.append((length,tip,engaged))
                if not options:
                    rows.append(dict(**row,hits=[dict(reason='no suitable bolt length')],preliminary_clear=False));continue
                length,tip,engaged=options[0];row.update(screw_length_mm=length,screw_tip_mm=tip,nominal_insert_engagement_mm=engaged)
                bolt=M.m3_parts(p,n,thickness,-60.,length)
                insert=rod(p+n*bottom,p+n*top,2.5)-rod(p+n*(bottom-.1),p+n*(top+.1),1.6)
                for name,body in [('insert',insert),('bolt',bolt['bolt']),('front_washer',bolt['front_washer'])]:
                    for key,obstacle in solids.items():
                        if key in ['bone:scapula','cap:FJ1467/origin']:continue
                        v=float(A.I.intersection_volume(body,obstacle))
                        if v<=.001:continue
                        if key=='tpu:deltoid':v=float(((body^obstacle)-zone).volume())
                        if v>.001:hits.append(dict(part=name,target=key,volume_mm3=v))
                    for head,wire in wires.items():
                        v=float(A.I.intersection_volume(body,wire))
                        if v>.001:hits.append(dict(part=name,target='active_wire:'+head,volume_mm3=v))
                rows.append(dict(**row,hits=hits,preliminary_clear=not hits))
        print('grid x',x,'eligible',sum(r['preliminary_clear'] for r in rows),flush=True)
    # Keep the shallowest feasible seat at a given axis; additional burial
    # does not justify a second axis at the same position.
    axes={}
    for row in rows:
        if row['preliminary_clear']:
            key=tuple(row['offset_local_mm'])
            if key not in axes or row['setback_mm']<axes[key]['setback_mm']:axes[key]=row
    pairs=[]
    for a,b in combinations(axes.values(),2):
        pa=np.array(a['offset_local_mm']);pb=np.array(b['offset_local_mm']);distance=float(np.linalg.norm(pa-pb))
        if distance<10.:continue
        pairs.append(dict(holes=[a,b],pitch_mm=distance,pair_midpoint_offset_mm=((pa+pb)*.5).tolist(),
            score=float(np.linalg.norm((pa+pb)*.5)+abs(distance-14.)*.25)))
    selected=min(pairs,key=lambda x:x['score']) if pairs else None
    save(DEST/'offset_insert_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        symmetric_study_sha256=sha(first_file),script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        grid_step_mm=2.,axis_radius_limit_mm=10.,insert_OD_mm=5.,insert_length_mm=6.,radial_stock_guard_diameter_mm=7.,
        minimum_axial_reserve_mm=1.,rows=rows,eligible_axes=len(axes),eligible_pairs=len(pairs),selected=selected,
        screen_pass=bool(selected),mount_geometry_generated=False,stage_pass=False,
        not_verified=['insert retention and real threads','actual bores and connected bone','TPU and support shapes','installation and actuation']))
    print('selected',selected,flush=True)


if __name__=='__main__':main()
