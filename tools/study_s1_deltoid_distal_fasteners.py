"""Screen shared deltoid distal through-fasteners in the current shoulder.

The three outlet centres remain fixed. Hole radius/rotation may change, not
the cap centre, bone placement or wire paths. Real cap/support construction
and assembly remain separate requirements after this envelope study.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_deltoid_distal_study'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_middle_side_nut_mount/assembly_check.json';parent=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['static_fit_stock_contact_pass'] and proof['local_approaches_geometry_pass']
    solids={k:load_component(s) for k,s in parent['sources'].items()};bone=solids['bone:humerus']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    port=next(r for r in routes if r['head']=='FJ1468')['ports']['distal']
    I=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-np.array(port['inner_face_source_mm']))@n)
    zone=rod(I-n*18.,I+n*.2,14.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    wall_guards={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),2.) for r in routes if r['head'] in ['FJ1468','FJ1467','FJ1513']}
    rows=[]
    for radius in [7.,8.,9.]:
        for angle in range(0,180,5):
            phi=np.deg2rad(angle);u=F[:,0]*np.cos(phi)+F[:,1]*np.sin(phi);holes=[];hits=[]
            for index,sign in enumerate([-1,1]):
                p=I+u*sign*radius;ray=bone.ray_cast(p+n*10.,p-n*65.)
                zs=[float((np.array(h.position)-I)@n) for h in ray]
                if len(zs)!=2:
                    hits.append(dict(index=index,reason='target is not a single bone segment',crossings=zs));continue
                guard=rod(p+n*(zs[1]+1.5),p+n*(zs[0]-1.5),2.7)
                missing=float((guard-bone).volume())
                if missing>.001:hits.append(dict(index=index,reason='bone radial stock guard not contained',missing_mm3=missing))
                holes.append(dict(index=index,axis_source_mm=p.tolist(),front_normal_mm=zs[0],rear_normal_mm=zs[1],
                    ray_stock_mm=zs[0]-zs[1],radial_stock_guard_missing_mm3=missing))
            if len(holes)!=2:
                rows.append(dict(radius_mm=radius,angle_deg=angle,holes=holes,hits=hits,preliminary_clear=False));continue
            lower=min(h['rear_normal_mm'] for h in holes)-2.2
            length=next(x for x in [35.,40.,45.,50.,55.,60.] if x>=thickness+.75-(lower-3.15)+.5)
            for index,h in enumerate(holes):
                p=np.array(h['axis_source_mm']);parts=M.m3_parts(p,n,thickness,lower,length)
                drill=rod(p-n*60.,p+n*.2,1.7)^zone
                for head,guard in wall_guards.items():
                    v=float(A.I.intersection_volume(drill,guard))
                    if v>.001:hits.append(dict(index=index,reason='new bore cuts protected wire wall',head=head,volume_mm3=v))
                for name,body in parts.items():
                    for key,obstacle in solids.items():
                        if key in ['bone:humerus','cap:FJ1468/distal']:continue
                        v=float(A.I.intersection_volume(body,obstacle))
                        if v<=.001:continue
                        if key=='tpu:deltoid':v=float(((body^obstacle)-zone).volume())
                        if v>.001:hits.append(dict(index=index,part=name,target=key,volume_mm3=v))
                    for head,wire in wires.items():
                        v=float(A.I.intersection_volume(body,wire))
                        if v>.001:hits.append(dict(index=index,part=name,target='active_wire:'+head,volume_mm3=v))
            rows.append(dict(radius_mm=radius,angle_deg=angle,holes=holes,rear_plate_lower_mm=lower,
                screw_length_mm=length,hits=hits,preliminary_clear=not hits))
        print('radius',radius,'clear',sum(r['preliminary_clear'] for r in rows if r['radius_mm']==radius),flush=True)
    eligible=[r for r in rows if r['preliminary_clear']]
    selected=min(eligible,key=lambda r:(r['radius_mm'],-min(h['ray_stock_mm'] for h in r['holes']))) if eligible else None
    save(DEST/'fastener_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        parent_verification_sha256=sha(proof_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        angle_step_deg=5.,hole_radii_mm=[7.,8.,9.],bone_stock_guard_radius_mm=2.7,
        permitted_terminal_radius_mm=14.,permitted_terminal_normal_mm=[-18.,.2],protected_three_wire_radius_mm=2.,
        rows=rows,selected=selected,screen_pass=bool(selected),mount_geometry_generated=False,stage_pass=False,
        not_verified=['actual holes and bone continuity','front/rear supports and TPU clamp geometry',
                      'local/whole assembly and tools','actual threading, preload and strength','actuation']))
    print('selected',None if not selected else (selected['radius_mm'],selected['angle_deg']),flush=True)


if __name__=='__main__':main()
