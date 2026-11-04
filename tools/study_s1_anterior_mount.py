"""Screen complete anterior M3 envelopes in the latest shoulder assembly.

Only a declared terminal zone may receive new TPU bores. Target clavicle and
original cap still require real bores/support geometry in a subsequent build.
No successful screen is a mount, assembly, or strength acceptance.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_anterior_mount_study'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_posterior_split_clamp/assembly_check.json'
    parent=read(file);proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['local_rigid_approaches_and_seats_pass']
    axis_file=A.DESIGN/'s1_deltoid_origin_holes/hole_axis_study.json';axis=read(axis_file)
    group=next(g for g in axis['groups'] if g['head']=='FJ1468')
    original_axis_bone=read(ROOT/axis['parent_report'])['sources']['bone:clavicle']
    assert original_axis_bone==parent['sources']['bone:clavicle']
    solids={k:load_component(s) for k,s in parent['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    port=next(r for r in routes if r['head']=='FJ1468')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source'])
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    zone=rod(I-n*18.,I+n*.2,12.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    rows=[]
    for pattern in group['rows']:
        if not pattern['both_axes_cross_target_once']:continue
        lower=min(h['bones']['bone:clavicle'][1]['normal_mm'] for h in pattern['holes'])-2.2
        for direction in ['front_head','rear_head']:
            if direction=='front_head':
                length=next(x for x in [35.,40.,45.,50.,55.,60.] if x>=thickness+.75-(lower-3.15)+.5)
            else:
                length=next(x for x in [35.,40.,45.,50.,55.,60.] if x>=2.+2.4-(lower-.75)+.5)
            hits=[]
            for index,hole in enumerate(pattern['holes']):
                point=np.array(hole['axis_source_mm'])
                if direction=='front_head':
                    parts=M.m3_parts(point,n,thickness,lower,length)
                else:
                    # Reuse actual M3 envelopes in reversed frame. Front nut
                    # starts at n=2mm inside a future keyed cap pocket.
                    parts=M.m3_parts(point,-n,-lower,-1.25,length)
                for name,body in parts.items():
                    for key,obstacle in solids.items():
                        if key in ['bone:clavicle','cap:FJ1468/origin']:continue
                        collision=body^obstacle
                        # A bore inside the bounded end boss is a planned
                        # change, not permission to drill working lattice.
                        if key=='tpu:deltoid':collision=collision-zone
                        volume=float(collision.volume())
                        if volume>.001:hits.append(dict(index=index,part=name,target=key,volume_mm3=volume))
                    for head,wire in wires.items():
                        volume=float(A.I.intersection_volume(body,wire))
                        if volume>.001:hits.append(dict(index=index,part=name,target='active_wire:'+head,volume_mm3=volume))
            rows.append(dict(angle_deg=pattern['angle_deg'],direction=direction,nominal_M3_length_mm=length,
                minimum_ray_stock_mm=pattern['minimum_stock_mm'],rear_plate_lower_mm=lower,
                holes=pattern['holes'],hits=hits,preliminary_clear=not hits))
            print(pattern['angle_deg'],direction,'hits',len(hits),flush=True)
    eligible=[r for r in rows if r['preliminary_clear']]
    selected=max(eligible,key=lambda r:r['minimum_ray_stock_mm']) if eligible else None
    save(DEST/'fastener_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),axis_study_sha256=sha(axis_file),
        routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),angle_step_deg=5.,
        cap_centre_source_mm=I.tolist(),normal_source=n.tolist(),cap_thickness_mm=thickness,
        permitted_terminal_zone=dict(radius_mm=12.,normal_interval_mm=[-18.,.2]),rows=rows,
        selected=selected,screen_pass=bool(eligible),mount_geometry_generated=False,stage_pass=False,
        not_verified=['cap and bone bores','front support and rear plate','TPU load-transfer flange',
                      'assembly and tool paths','threads, preload and strength','driven range']))
    print('selected',None if selected is None else (selected['angle_deg'],selected['direction']),flush=True)


if __name__=='__main__':main()
