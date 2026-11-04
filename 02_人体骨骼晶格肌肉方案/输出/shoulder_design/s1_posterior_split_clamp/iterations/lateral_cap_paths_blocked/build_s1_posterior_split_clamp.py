"""Posterior split cap with front antirotation nut pockets and rear-fed screws.

Nuts are manually held in open pockets during local assembly; permanent
capture, preload, material strength and actual tool engagement are unverified.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,sha,read,save,export,load_component,shoulder_pads

DEST=A.DESIGN/'s1_posterior_split_clamp'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_deltoid_posterior_mount/assembly_check.json';parent=read(file)
    proof_file=file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['static_fit_and_contact_pass']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    route=next(row for row in routes if row['head']=='FJ1513');port=route['ports']['origin']
    I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);P=shoulder_pads.pad_basis(n)
    direction=np.array(route['core_path_source_mm'][0])-I
    tangent=A.D.unit(direction-n*(direction@n));side=A.D.unit(np.cross(tangent,n))
    radians=np.deg2rad(parent['dimensions']['hole_rotation_deg'])
    hole_axis=P[:,0]*np.cos(radians)+P[:,1]*np.sin(radians)
    if abs(side@hole_axis)<.75:
        side=hole_axis.copy();tangent=np.cross(n,side)
        if tangent@direction<0:side=-side;tangent=-tangent
    F=np.stack([side,tangent,n],axis=1);mf=A.I._manifold()
    old_cap=load_component(parent['sources']['cap:FJ1513/origin'])
    cap=old_cap;parts={};holes=[]
    nut_base=2.;nut_height=2.4
    # The original central cap disc is smaller than its nominal 14mm hole
    # pattern. Add real nut-bearing lugs; do not place pockets in empty space.
    for hole in parent['selected_bone_crossings']:
        point=np.array(hole['axis_source_mm'])
        cap=cap+rod(point,point+n*4.8,3.7)
        cap=cap-rod(point-n*.2,point+n*5.,1.7)
    for index,hole in enumerate(parent['selected_bone_crossings']):
        point=np.array(hole['axis_source_mm']);holes.append(point)
        pocket=mf.Manifold.cylinder(3.,5.7/np.sqrt(3),5.7/np.sqrt(3),6).translate([0,0,nut_base])
        pocket=pocket.transform(np.c_[P,point]);cap=cap-pocket
        nut=mf.Manifold.cylinder(nut_height,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,0,nut_base])
        nut=nut-mf.Manifold.cylinder(2.8,1.6,1.6,64).translate([0,0,nut_base-.2])
        parts[f'origin:FJ1513/nut_{index}']=nut.transform(np.c_[P,point])
        head_base=parent['dimensions']['rear_plate_lower_mm']-.75
        length=25.
        profile=np.array([[0,-length],[1.5,-length],[1.5,0],[2.75,0],[2.75,3.],[0,3.]])
        bolt=mf.Manifold.revolve(mf.CrossSection([profile]),64)
        socket=mf.Manifold.cylinder(1.7,2.5/np.sqrt(3),2.5/np.sqrt(3),6).translate([0,0,1.5])
        bolt=bolt-socket
        rear_frame=shoulder_pads.pad_basis(-n)
        parts[f'origin:FJ1513/bolt_{index}']=bolt.transform(np.c_[rear_frame,point+n*head_base])
    local=cap.transform(np.c_[F.T,-F.T@I])
    left=local^mf.Manifold.cube([19.9,40.,20.]).translate([-20.,-20.,-5.])
    right=local^mf.Manifold.cube([19.9,40.,20.]).translate([.1,-20.,-5.])
    for y in [-2.5,-4.5]:
        left=left+rod([-1.1,y,1.5],[1.2,y,1.5],.65)
        right=right-rod([-1.4,y,1.5],[1.4,y,1.5],.75)
    parts['cap:FJ1513/origin_A']=left.transform(np.c_[F,I])
    parts['cap:FJ1513/origin_B']=right.transform(np.c_[F,I])
    original={key:load_component(row) for key,row in parent['sources'].items()}
    contour=[]
    for key in ['cap:FJ1513/origin_A','cap:FJ1513/origin_B']:
        before=parts[key]
        for target in ['tpu:deltoid','tpu:infraspinatus']:
            obstacle=original[target];guard=obstacle
            for axis in range(3):
                for sign in [-1,1]:guard=guard+obstacle.translate(np.eye(3)[axis]*sign*.03)
            parts[key]=parts[key]-guard
        contour.append(dict(key=key,removed_mm3=float((before-parts[key]).volume())))
    retired=['cap:FJ1513/origin','origin:FJ1513/front_washer_0','origin:FJ1513/front_washer_1']
    replaced=[f'origin:FJ1513/{kind}_{index}' for kind in ['bolt','nut'] for index in range(2)]
    sources={key:row for key,row in parent['sources'].items() if key not in retired+replaced}
    solids={key:original[key] for key in sources}
    cleanup=[]
    for key,solid in parts.items():
        prepared=solid.as_original().simplify(1e-6)
        cleanup.append(dict(key=key,added_mm3=float((prepared-solid).volume()),removed_mm3=float((solid-prepared).volume())))
        solids[key],sources[key]=export(prepared,DEST/(key.replace(':','_').replace('/','_')+'.obj'),I)
    hits=[];count=0;keys=list(solids)
    for index,a in enumerate(keys):
        for b in keys[index+1:]:
            if a not in parts and b not in parts:continue
            count+=1;v=float(A.I.intersection_volume(solids[a],solids[b]))
            old=float(A.I.intersection_volume(original[a],original[b])) if a in original and b in original else 0.
            if v>max(.001,old+.001):hits.append(dict(a=a,b=b,volume_mm3=v,parent_mm3=old))
    wire_checks=[]
    for row in routes:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        wh=[]
        for key in parts:
            v=float(A.I.intersection_volume(solids[key],wire))
            if v>.001:wh.append(dict(key=key,volume_mm3=v))
        wire_checks.append(dict(head=row['head'],hits=wh))
    assignments={}
    for index,point in enumerate(holes):assignments[str(index)]='A' if (point-I)@side<0 else 'B'
    assert set(assignments.values())=={'A','B'}
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(file),parent_verification_sha256=sha(proof_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        retired_keys=retired,replaced_keys=replaced,new_keys=['cap:FJ1513/origin_A','cap:FJ1513/origin_B'],
        local_origin_source_mm=I.tolist(),local_basis_columns=F.tolist(),nut_to_half=assignments,
        dimensions=dict(split_gap_mm=.2,locator_diameter_mm=1.3,locator_pocket_diameter_mm=1.5,
            locator_tangent_mm=[-2.5,-4.5],locator_axial_mm=1.5,nut_pocket_AF_mm=5.7,
            nut_lug_radius_mm=3.7,nut_lug_height_mm=4.8,
            nut_floor_mm=2.,nut_height_mm=2.4,nut_nominal_AF_mm=5.5,
            rear_M3_length_mm=25.,rear_head_base_normal_mm=head_base,bolt_front_tip_normal_mm=head_base+25.),
        unchanged_parent_rows=len(parent['sources'])-len(retired)-len(replaced),numerical_cleanup=cleanup,
        muscle_contour_relief=contour,contour_guard_axis_shifts_mm=.03,
        affected_pairs=count,static_hits=hits,six_active_wire_corridors=wire_checks,
        static_geometry_pass=not hits and not any(row['hits'] for row in wire_checks),
        manual_nut_retention_during_assembly_required=True,permanent_nut_capture_verified=False,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**{k:v for k,v in parent['intended_body_owners'].items() if k not in retired},
            'cap:FJ1513/origin_A':'scapula_ideal_fastened','cap:FJ1513/origin_B':'scapula_ideal_fastened'},
        adopted_as_current=False,stage_pass=False,whole_module_installation_verified=False,actual_drive_range_deg=None)
    save(DEST/'assembly_check.json',report)
    print(dict(parts=len(sources),pairs=count,hits=hits,wires=[row for row in wire_checks if row['hits']]),flush=True)
    if not report['static_geometry_pass']:raise ValueError('Posterior split clamp needs revision; saved evidence retained')


if __name__=='__main__':main()
