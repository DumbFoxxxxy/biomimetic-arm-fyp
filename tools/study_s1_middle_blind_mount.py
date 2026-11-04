"""Screen a front-installed blind-insert middle-deltoid connection.

Multiple deeper scapular ray crossings do not invalidate a blind bore: only
the first contiguous surface segment is used, with a full radial stock guard.
Nominal insert OD5 x 6 and M3 envelopes are design variables, not a sourced
part, thread, retention, or strength validation. No bone is modified here.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component,shoulder_pads

DEST=A.DESIGN/'s1_middle_blind_study'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_anterior_mount/assembly_check.json';parent=read(file)
    proof_file=file.parent/'side_entry_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['static_fit_and_contact_pass'] and proof['local_approaches_geometry_pass']
    solids={k:load_component(s) for k,s in parent['sources'].items()}
    bone=solids['bone:scapula'];mf=A.I._manifold()
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    port=next(r for r in routes if r['head']=='FJ1467')['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
    thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    zone=rod(I-n*18.,I+n*.2,12.)
    wires={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3) for r in routes}
    rows=[]
    for angle in range(0,180,5):
        phi=np.deg2rad(angle);u=F[:,0]*np.cos(phi)+F[:,1]*np.sin(phi)
        for setback in [1.,2.,3.]:
            holes=[];hits=[]
            for index,sign in enumerate([-1,1]):
                point=I+u*sign*7.;crossings=bone.ray_cast(point+n*10.,point-n*80.)
                if len(crossings)<2:
                    hits.append(dict(index=index,reason='no first bone segment'));continue
                zs=[float((np.array(h.position)-I)@n) for h in crossings]
                top=zs[0]-setback;bottom=top-6.
                guard=rod(point+n*bottom,point+n*top,3.5)
                missing=float((guard-bone).volume())
                if missing>.001:
                    hits.append(dict(index=index,reason='radial/bottom guard not contained',missing_mm3=missing))
                if bottom-.5<zs[1]+1.:
                    hits.append(dict(index=index,reason='less than 1mm axial reserve',reserve_mm=bottom-.5-zs[1]))
                options=[]
                for length in [10.,12.,14.,16.,18.,20.]:
                    tip=thickness+.75-length;engaged=top-max(bottom,tip)
                    if engaged>=4. and tip>=bottom-.3:
                        options.append((length,tip,engaged))
                if not options:
                    hits.append(dict(index=index,reason='no length with 4mm insert engagement'));continue
                length,tip,engaged=options[0]
                bolt_parts=M.m3_parts(point,n,thickness,-60.,length)
                insert=rod(point+n*bottom,point+n*top,2.5)-rod(point+n*(bottom-.1),point+n*(top+.1),1.6)
                for name,body in [('insert',insert),('bolt',bolt_parts['bolt']),('front_washer',bolt_parts['front_washer'])]:
                    for key,obstacle in solids.items():
                        if key in ['bone:scapula','cap:FJ1467/origin']:continue
                        v=float(A.I.intersection_volume(body,obstacle))
                        if v<=.001:continue
                        if key=='tpu:deltoid':v=float(((body^obstacle)-zone).volume())
                        if v>.001:hits.append(dict(index=index,part=name,target=key,volume_mm3=v))
                    for head,wire in wires.items():
                        v=float(A.I.intersection_volume(body,wire))
                        if v>.001:hits.append(dict(index=index,part=name,target='active_wire:'+head,volume_mm3=v))
                holes.append(dict(index=index,axis_source_mm=point.tolist(),ray_crossings_normal_mm=zs,
                    first_segment_stock_mm=zs[0]-zs[1],insert_top_mm=top,insert_bottom_mm=bottom,
                    radial_guard_missing_mm3=missing,bottom_axial_reserve_mm=bottom-.5-zs[1],
                    screw_length_mm=length,screw_tip_mm=tip,nominal_insert_engagement_mm=engaged))
            rows.append(dict(angle_deg=angle,setback_mm=setback,holes=holes,hits=hits,
                preliminary_clear=len(holes)==2 and not hits))
        print(angle,'clear',sum(r['preliminary_clear'] for r in rows if r['angle_deg']==angle),flush=True)
    eligible=[r for r in rows if r['preliminary_clear']]
    selected=max(eligible,key=lambda r:min(h['bottom_axial_reserve_mm'] for h in r['holes'])) if eligible else None
    save(DEST/'blind_insert_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        angle_step_deg=5.,pitch_mm=14.,insert_outer_diameter_mm=5.,insert_length_mm=6.,
        radial_guard_diameter_mm=7.,minimum_axial_reserve_mm=1.,rows=rows,selected=selected,
        screen_pass=bool(eligible),mount_geometry_generated=False,stage_pass=False,
        not_verified=['blind bore and solid continuity after cutting','insert retention method','support and TPU clamp geometry',
                      'local and whole installation paths','physical preload and pullout strength','actuation']))
    print('selected',None if not selected else (selected['angle_deg'],selected['setback_mm']),flush=True)


if __name__=='__main__':main()
