"""Front supports/backing plates for the two proximal shoulder muscle caps.

Actual drilled blade, four nominal M3 screws, nuts and washers; TPU clamp
pads remain part of the muscles. All available parts are checked, with no
whole-pair bone/contact exclusions. No strength/preload/drive claim.
"""
import itertools
import json
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import design_shoulder_tendon_link as L
import shoulder_pads
from build_shoulder_shaft_collar import rod
from build_shoulder_isolated_retained_joint import export_checked
from build_shoulder_unified_muscle import local_export

SOURCE=L.DEST/'assembly_check.json'
DEST=L.BASE/'origin_mounts'


def m3_parts(inner,normal,cap_thickness,back_outer,length):
    """Nominal plain envelopes, no thread engagement or actual preload."""
    mf=A.I._manifold();n=A.D.unit(normal);F=shoulder_pads.pad_basis(n)
    head=cap_thickness+.75
    profile=np.array([[0,head-length],[1.5,head-length],[1.5,head],
                      [2.75,head],[2.75,head+3],[0,head+3]])
    bolt=mf.Manifold.revolve(mf.CrossSection([profile]),64).transform(np.c_[F,inner])
    nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,0,back_outer-3.15])
    nut=nut-mf.Manifold.cylinder(4.,1.6,1.6,48).translate([0,0,back_outer-4])
    parts=dict(bolt=bolt,nut=nut.transform(np.c_[F,inner]))
    for key,a,b in [('front_washer',cap_thickness+.15,cap_thickness+.65),
                    ('back_washer',back_outer-.65,back_outer-.15)]:
        parts[key]=(rod(inner+n*a,inner+n*b,3.5)-
                    rod(inner+n*(a-.2),inner+n*(b+.2),1.7))
    return parts


def clearance(solid,amount=.1):
    parts=[solid]
    parts.extend(solid.translate(np.eye(3)[axis]*sign*amount)
                 for axis in range(3) for sign in [-1,1])
    return A.I._manifold().Manifold.batch_boolean(parts,A.I._manifold().OpType.Add)


def remove_obsolete_terminal_residue(solid,zone,channel):
    """Only isolated remnants of the replaced terminal, never working lattice.

Require the explicit terminal zone, no channel-wall intersection, <1mm3 total
and <0.01% of the body. A larger/differently located fragment stops the build.
The existing cavity shells are retained by Boolean difference, not compose.
"""
    positive=sorted([p for p in solid.decompose() if p.volume()>1e-9],key=lambda p:-p.volume())
    if len(positive)<=1:return solid,[]
    tiny=positive[1:];total=sum(p.volume() for p in tiny)
    if total>=1. or total/max(1.,solid.volume())>=1e-4:
        raise ValueError('New mount interrupts useful TPU; redesign instead of discarding pieces')
    removed=[]
    for fragment in tiny:
        outside=float((fragment-zone).volume())
        protected=A.I.intersection_volume(fragment,channel)
        if outside>1e-7 or protected>1e-7:
            raise ValueError('Fragment is outside the replaced terminal or intersects its channel')
        removed.append(dict(volume_mm3=float(fragment.volume()),bounds_source_mm=list(fragment.bounding_box()),
            outside_replaced_zone_mm3=outside,protected_channel_intersection_mm3=protected,
            reason='isolated obsolete soft terminal remnant after replacement by actual saddle/pad'))
        solid=solid-fragment
    return solid,removed


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    parent=json.loads(SOURCE.read_text(encoding='utf8'))
    audit_file=L.BASE/'origin_fixation_review/fixation_audit.json'
    audit=json.loads(audit_file.read_text(encoding='utf8'))
    if audit['source_sha256']!=A.hashed(SOURCE):raise ValueError('Stale proximal ray audit')
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes']
    solids={k:A.load_checked(A.ROOT/v['file'],v['translation_source_mm'],v['sha256'])[0]
            for k,v in parent['sources'].items()}
    original_bone=solids['bone:scapula'];drills=[];builds={};new={};new_meta={};replacement_tpu={}
    guards=[A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),1.95) for r in routes]
    for row in audit['rows']:
        head=row['head'];route=next(r for r in routes if r['head']==head)
        port=route['ports']['origin'];n=A.D.unit(port['normal_source']);I=np.array(port['inner_face_source_mm'])
        F=shoulder_pads.pad_basis(n);cap_thickness=float(np.dot(np.array(port['outer_face_source_mm'])-I,n))
        pad_height=2.5 if head=='FJ1500' else 4.;base_height=2.
        front=rod(I-n*(pad_height+base_height),I-n*pad_height,11.)
        back_face=min(h['bone_crossings'][1]['position_along_cap_normal_mm'] for h in row['holes'])-.2
        back_outer=back_face-2.;parts={};centres=[]
        for hole in row['holes']:
            sign=hole['sign'];Q=np.array(hole['axis_point_source_mm']);centres.append(Q)
            z0=hole['bone_crossings'][0]['position_along_cap_normal_mm']
            z1=hole['bone_crossings'][1]['position_along_cap_normal_mm']
            if hole['bone_stock_mm']<=1.1:raise ValueError('Foot construction needs review for thinner blade')
            front=front+rod(Q+n*(z0-1),Q+n*(-pad_height-.1),3.8)
            drill=rod(Q+n*10,Q-n*60,1.7);drills.append(drill)
        # A narrow capsule bridges two bearing feet; it is not a full R11
        # disc across unrelated muscle origins on the opposite blade surface.
        back=(rod(centres[0]+n*back_outer,centres[0]+n*back_face,4.)+
              rod(centres[1]+n*back_outer,centres[1]+n*back_face,4.))
        bridge=A.I._manifold().Manifold.cube([14.,8.,2.]).translate([-7.,-4.,back_outer])
        back=back+bridge.transform(np.c_[F,I])
        for hole,Q in zip(row['holes'],centres):
            z1=hole['bone_crossings'][1]['position_along_cap_normal_mm']
            back=back+rod(Q+n*(back_face-.1),Q+n*(z1+1),3.8)
        front=front-original_bone;back=back-original_bone
        for Q in centres:
            bore=rod(Q+n*10,Q-n*60,1.7);front=front-bore;back=back-bore
        for guard in guards:
            if A.I.intersection_volume(front,guard)>1e-9:front=front-guard
            if A.I.intersection_volume(back,guard)>1e-9:back=back-guard
        parts.update(front_support=front,backing_plate=back)
        nut_start=back_outer-3.15
        minimum=cap_thickness+.75-nut_start+1.
        length=float(next(v for v in [16,20,25,30,35,40,45] if v>=minimum))
        for index,Q in enumerate(centres):
            parts.update({key+'_'+str(index):v for key,v in m3_parts(Q,n,cap_thickness,back_outer,length).items()})
        muscle_key='tpu:'+route['module'];old_tpu=solids[muscle_key]
        pad=rod(I-n*pad_height,I,11.)
        # Keep the working belly and its channel, add the actual flat clamp
        # pad and relieve the real saddle. Only bounded obsolete terminal
        # remnants can be removed, under the explicit checks below.
        zone=rod(I-n*8.,I+n*.01,11.2)
        lumen=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.86)
        wall=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),1.8)
        # Intentional flat clamp contact is not a fit gap: relieve the old
        # belly around the support, then place the new pad directly on its
        # upper face. Post ends stop inside the flange, not in the soft pad.
        revised=(old_tpu-clearance(front))+pad
        for Q in centres:
            revised=revised-rod(Q-n*(pad_height+3),Q+n*(cap_thickness+2),1.7)
        revised=revised-lumen
        revised,residue=remove_obsolete_terminal_residue(revised,zone,wall)
        folder=DEST/head;folder.mkdir(exist_ok=True)
        if (folder/'rejected_TPU.json').exists():
            archive=folder/'iterations/initial_pad_interrupted';archive.mkdir(parents=True,exist_ok=True)
            for filename in ['rejected_TPU.json','TPU_proximal_clamp_pad.stl']:
                if (folder/filename).exists() and not (archive/filename).exists():
                    shutil.copyfile(folder/filename,archive/filename)
        original=parent['sources'][muscle_key]
        actual,meta=local_export(revised,folder/'TPU_proximal_clamp_pad.stl',np.array(original['translation_source_mm']))
        if not meta['mesh_checks_pass']:
            (folder/'rejected_TPU.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf8',newline='\n')
            raise ValueError('Proximal TPU readback/connectedness fails '+head)
        solids[muscle_key]=actual;replacement_tpu[muscle_key]=meta
        for key,body in parts.items():
            label='origin:'+head+'/'+key
            if key in ['front_support','backing_plate']:
                # Keep these small conforming parts near their own origin;
                # the anatomical display origin needlessly loses float32
                # precision at the thin foot/bone surface intersections.
                actual,meta=local_export(body.as_original().simplify(.001),folder/(key+'.stl'),I)
                if not meta['mesh_checks_pass']:raise ValueError('Support STL readback fails '+head+'/'+key)
            else:
                actual,meta=export_checked(body.as_original().simplify(.001),folder/(key+'.stl'))
            new[label]=actual;new_meta[label]=meta
        builds[head]=dict(pad_height_mm=pad_height,back_outer_from_cap_inner_mm=back_outer,
            nominal_M3_length_mm=length,front_support_volume_mm3=float(front.volume()),
            backing_volume_mm3=float(back.volume()),TPU_added_mm3=float((solids[muscle_key]-old_tpu).volume()),
            TPU_removed_mm3=float((old_tpu-solids[muscle_key]).volume()),
            terminal_cleanup_zone=dict(radius_mm=11.2,normal_interval_from_cap_inner_mm=[-8.,.01]),
            removed_obsolete_terminal_residue=residue,
            additional_missing_protected_channel_wall_mm3=float((
                ((A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),1.61)-lumen)^old_tpu)-solids[muscle_key]).volume()),
            nominal_plate_contact_not_preload_test=True)
        print(head,'M3x',length,'front/back',round(front.volume(),2),round(back.volume(),2),flush=True)
    drill_union=A.I._manifold().Manifold.batch_boolean(drills,A.I._manifold().OpType.Add)
    drilled,drilled_meta=export_checked((original_bone-drill_union).as_original().simplify(.001),DEST/'scapula_four_M3_holes.stl')
    solids['bone:scapula']=drilled
    sources=dict(parent['sources']);sources['bone:scapula']=drilled_meta;sources.update(replacement_tpu);sources.update(new_meta)
    solids.update(new);contacts=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:contacts.append(dict(a=a,b=b,volume_mm3=volume,
            category='compliant contact review' if a.startswith('tpu:') or b.startswith('tpu:')
            or a.startswith('collar:liner') or b.startswith('collar:liner') else 'rigid/wire fit review'))
    wires=[]
    for row in routes:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(wire,s)) for k,s in solids.items()]
        wires.append(dict(head=row['head'],hits=hits,hits_over_threshold=[h for h in hits if h['volume_mm3']>.001]))
    owners=dict(parent['intended_body_owners']);owners.update({k:'scapula_ideal_fastened' for k in new_meta})
    report=dict(scope=__doc__,n_instances=len(sources),sources=sources,
        parent_source=SOURCE.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(SOURCE),
        audit_sha256=A.hashed(audit_file),routes_sha256=A.hashed(route_file),script_sha256=A.hashed(Path(__file__).resolve()),
        proximal_mounts=builds,bone_drilling=dict(holes=4,radius_mm=1.7,
            removed_from_parent_mm3=float((original_bone-drilled).volume()),added_by_export_mm3=float((drilled-original_bone).volume())),
        material_intersections=contacts,actuation_wire_checks=wires,intended_body_owners=owners,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        not_verified=['M3 torque/preload, blade bearing, creep, friction and strength',
            'assembly tool access and actual pad compression',
            'same unresolved flexible muscle/tendon equilibrium and drive range as parent',
            'no physical printing or axis-coupon measurements'],stage_pass=False,actual_drive_range_deg=None)
    (DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Available objects',len(sources),'rigid hits',
        [r for r in contacts if r['volume_mm3']>.001 and r['category'].startswith('rigid')],
        'wire hits',sum(bool(w['hits_over_threshold']) for w in wires),flush=True)
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
