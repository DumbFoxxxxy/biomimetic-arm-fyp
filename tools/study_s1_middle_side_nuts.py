"""Middle deltoid: shallow side-loaded hex nuts instead of blind inserts.

Roof and floor stock are checked in full solid envelopes. Later a builder
must cut and verify real side slots and insertion paths. Nominal threads and
stock checks are not pullout strength or preload verification.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_middle_blind_study'


def main():
    offset_file=DEST/'offset_insert_study.json';offset=read(offset_file);assert not offset['screen_pass']
    file=ROOT/offset['parent_report'];parent=read(file);assert sha(file)==offset['parent_sha256']
    solids={k:load_component(s) for k,s in parent['sources'].items()};bone=solids['bone:scapula'];mf=A.I._manifold()
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    port=next(r for r in routes if r['head']=='FJ1467')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n);zone=rod(I-n*18.,I+n*.2,12.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    rows=[];axes={}
    for x in range(-10,11,2):
        for y in range(-10,11,2):
            if x*x+y*y>100:continue
            p=I+F[:,0]*x+F[:,1]*y;ray=bone.ray_cast(p+n*10.,p-n*80.)
            zs=[float((np.array(h.position)-I)@n) for h in ray]
            if len(zs)<2:continue
            for setback in [2.,3.,4.,5.]:
                top=zs[0]-setback;bottom=top-2.4
                # R4.3 exceeds the AF5.7 pocket circumradius by ~1mm.
                guard=rod(p+n*(bottom-1.),p+n*(top+1.5),4.3)
                missing=float((guard-bone).volume());hits=[]
                row=dict(offset_local_mm=[x,y],axis_source_mm=p.tolist(),setback_mm=setback,
                    nut_top_mm=top,nut_bottom_mm=bottom,stock_guard_missing_mm3=missing,ray_crossings_normal_mm=zs)
                if missing>.001:
                    rows.append(dict(**row,hits=[dict(reason='roof/floor/radial guard exits bone')],preliminary_clear=False));continue
                options=[]
                for length in [10.,12.,14.,16.,18.,20.,22.,25.]:
                    tip=thickness+.75-length
                    if tip<=bottom-.2 and tip>=zs[1]+1.5:options.append((length,tip))
                if not options:
                    rows.append(dict(**row,hits=[dict(reason='no full nut engagement with bottom reserve')],preliminary_clear=False));continue
                length,tip=options[0];parts=M.m3_parts(p,n,thickness,-60.,length)
                nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,0,bottom])
                nut=nut-mf.Manifold.cylinder(2.8,1.6,1.6,48).translate([0,0,bottom-.2])
                nut=nut.transform(np.c_[F,p]);row.update(screw_length_mm=length,screw_tip_mm=tip)
                for name,body in [('nut',nut),('bolt',parts['bolt']),('front_washer',parts['front_washer'])]:
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
                if not hits and (x,y) not in axes:axes[(x,y)]=rows[-1]
        print('grid x',x,'eligible axes',len(axes),flush=True)
    pairs=[]
    for a,b in combinations(axes.values(),2):
        pa=np.array(a['offset_local_mm']);pb=np.array(b['offset_local_mm']);distance=float(np.linalg.norm(pa-pb))
        if distance<9.:continue
        pairs.append(dict(holes=[a,b],pitch_mm=distance,pair_midpoint_offset_mm=((pa+pb)*.5).tolist(),
            score=float(np.linalg.norm((pa+pb)*.5)+abs(distance-14.)*.25)))
    selected=min(pairs,key=lambda x:x['score']) if pairs else None
    save(DEST/'side_nut_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        failed_insert_study_sha256=sha(offset_file),script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        grid_step_mm=2.,axis_radius_limit_mm=10.,nut_AF_mm=5.5,nut_height_mm=2.4,pocket_AF_mm=5.7,
        stock_guard_diameter_mm=8.6,roof_stock_mm=1.5,floor_stock_mm=1.,rows=rows,
        eligible_axes=len(axes),eligible_pairs=len(pairs),selected=selected,screen_pass=bool(selected),stage_pass=False,
        not_verified=['cut side slot and retained roof/floor','nut insertion and lateral escape stops','TPU/support/cap shape',
                      'assembly, preload, strength and actuation']))
    print('selected',selected,flush=True)


if __name__=='__main__':main()
