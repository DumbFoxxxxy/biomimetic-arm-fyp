"""Offset shared-deltoid fasteners around the three preserved wire ports.

Independent fastener positions may move within a declared enlarged terminal
footprint. No source geometry is modified; supports, real bores, protection
and assembly must subsequently be proved on exported geometry.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_deltoid_distal_study'


def main():
    first_file=DEST/'fastener_study.json';first=read(first_file);assert not first['screen_pass']
    file=ROOT/first['parent_report'];parent=read(file);assert sha(file)==first['parent_sha256']
    solids={k:load_component(s) for k,s in parent['sources'].items()};bone=solids['bone:humerus']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    assert sha(route_file)==first['routes_sha256']
    port=next(r for r in routes if r['head']=='FJ1468')['ports']['distal']
    I=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-np.array(port['inner_face_source_mm']))@n)
    zone=rod(I-n*18.,I+n*.2,18.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    guards={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),2.) for r in routes if r['head'] in ['FJ1468','FJ1467','FJ1513']}
    rows=[]
    for x in range(-14,15,2):
        for y in range(-10,11,2):
            if x*x+y*y>196:continue
            p=I+F[:,0]*x+F[:,1]*y;ray=bone.ray_cast(p+n*10.,p-n*65.);zs=[float((np.array(h.position)-I)@n) for h in ray]
            hits=[];row=dict(offset_local_mm=[x,y],axis_source_mm=p.tolist())
            if len(zs)!=2:
                rows.append(dict(**row,hits=[dict(reason='not a single bone segment')],preliminary_clear=False));continue
            guard=rod(p+n*(zs[1]+1.5),p+n*(zs[0]-1.5),2.7);missing=float((guard-bone).volume())
            row.update(front_normal_mm=zs[0],rear_normal_mm=zs[1],ray_stock_mm=zs[0]-zs[1],radial_stock_guard_missing_mm3=missing)
            if missing>.001:
                rows.append(dict(**row,hits=[dict(reason='bone radial guard not contained',missing_mm3=missing)],preliminary_clear=False));continue
            lower=zs[1]-2.2;length=next(t for t in [35.,40.,45.,50.,55.,60.] if t>=thickness+.75-(lower-3.15)+.5)
            row.update(rear_plate_lower_mm=lower,screw_length_mm=length)
            drill=rod(p-n*60.,p+n*.2,1.7)^zone
            for head,guard in guards.items():
                v=float(A.I.intersection_volume(drill,guard))
                if v>.001:hits.append(dict(reason='cuts R2 wire guard',head=head,volume_mm3=v))
            for name,body in M.m3_parts(p,n,thickness,lower,length).items():
                for key,obstacle in solids.items():
                    if key in ['bone:humerus','cap:FJ1468/distal']:continue
                    v=float(A.I.intersection_volume(body,obstacle))
                    if v<=.001:continue
                    if key=='tpu:deltoid':v=float(((body^obstacle)-zone).volume())
                    if v>.001:hits.append(dict(part=name,target=key,volume_mm3=v))
                for head,wire in wires.items():
                    v=float(A.I.intersection_volume(body,wire))
                    if v>.001:hits.append(dict(part=name,target='active_wire:'+head,volume_mm3=v))
            rows.append(dict(**row,hits=hits,preliminary_clear=not hits))
        print('grid x',x,'eligible',sum(r['preliminary_clear'] for r in rows),flush=True)
    eligible=[r for r in rows if r['preliminary_clear']];pairs=[]
    for a,b in combinations(eligible,2):
        pa=np.array(a['offset_local_mm']);pb=np.array(b['offset_local_mm']);pitch=float(np.linalg.norm(pa-pb))
        if pitch<8.:continue
        # A bridge cap can use offset holes. Prefer the least footprint
        # extension, then retain a useful screw spacing near the old 14mm.
        score=float(max(np.linalg.norm(pa),np.linalg.norm(pb))+abs(pitch-14.)*.15+np.linalg.norm((pa+pb)*.5)*.1)
        pairs.append(dict(holes=[a,b],pitch_mm=pitch,score=score))
    selected=min(pairs,key=lambda x:x['score']) if pairs else None
    save(DEST/'offset_fastener_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        symmetric_study_sha256=sha(first_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        grid_step_mm=2.,axis_radius_limit_mm=14.,terminal_zone_radius_mm=18.,terminal_normal_mm=[-18.,.2],
        bone_guard_radius_mm=2.7,protected_three_wire_radius_mm=2.,rows=rows,eligible_axes=len(eligible),
        eligible_pairs=len(pairs),selected=selected,screen_pass=bool(selected),mount_geometry_generated=False,stage_pass=False))
    print('selected',selected,flush=True)


if __name__=='__main__':main()
