"""Teres active-wire hexagonal metal stop, fixed-cap bearing candidate.

External hex permits counterholding. Smooth M2x3 first-touch envelope does
not establish friction, actual thread engagement, torque or steel grip.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import (
    ROOT,A,M,shoulder_pads,rod,sha,read,save,export,load_component,
)
import design_shoulder_tendon_link as L

DEST=A.DESIGN/'s1_teres_wire_stop'


def main():
    DEST.mkdir(exist_ok=True)
    parent_file=A.DESIGN/'s1_teres_collar_adapter/assembly_check.json'
    parent=read(parent_file);proof_file=parent_file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(parent_file) and proof['local_interface_geometry_pass'] and proof['bare_collar_closing_path_pass']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(route_file)['routes'];route=next(row for row in routes if row['head']=='FJ1507')
    I=np.array(route['ports']['distal']['inner_face_source_mm']);n=A.D.unit(route['ports']['distal']['normal_source'])
    pad=shoulder_pads.pad_basis(n)
    # Keep the radial driver away from the M3 pair on pad-local X.
    x=pad[:,1];y=np.cross(n,x);F=np.stack([x,y,n],axis=1)
    affine=np.c_[F,I];mf=A.I._manifold()
    cylinder=lambda a,b,r:mf.Manifold.cylinder(b-a,r,r,96).translate([0,0,a])
    # Bearing head seats on the actual OUTER face; no source counterbore is
    # assumed. A short narrow pilot locates in the original wire bore.
    head=mf.Manifold.cylinder(4.,5./np.sqrt(3),5./np.sqrt(3),6).translate([0,0,4.05])
    blank=head+cylinder(3.05,4.05,.7)
    through=cylinder(2.8,8.3,.4)
    radial_F=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
    radial=lambda solid:solid.transform(np.c_[radial_F,[0.,0.,6.05]])
    port=radial(mf.Manifold.cylinder(4.,1.05,1.05,64))
    local_body=blank-through-port
    screw=(cylinder(.3,.6,1.) ^ mf.Manifold.cylinder(.3,.45,1.,64).translate([0,0,.3]))+cylinder(.6,3.3,1.)
    socket=mf.Manifold.cylinder(1.,.9/np.sqrt(3),.9/np.sqrt(3),6).translate([0,0,2.5])
    local_screw=radial(screw-socket)
    path=np.array(route['full_wire_path_source_mm']);tail=I+n*9.05
    extended=np.vstack([path,tail])
    objects={'stop:FJ1507/hex_body':local_body.transform(affine),
             'stop:FJ1507/M2x3':local_screw.transform(affine),
             'wire_active:FJ1507/reference':L.wire_tube(extended)}
    sources=dict(parent['sources']);solids={k:load_component(v) for k,v in sources.items()}
    for key,solid in objects.items():
        name=key.split('/')[-1]
        solids[key],sources[key]=export(solid,DEST/(name+'.obj'),I)
    pairs=[];count=0;keys=list(solids)
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            if a not in objects and b not in objects:continue
            count+=1;v=float(A.I.intersection_volume(solids[a],solids[b]))
            if v>.001:pairs.append(dict(a=a,b=b,volume_mm3=v))
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=parent_file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(parent_file),parent_verification_sha256=sha(proof_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        new_keys=sorted(objects),replaced_keys=[],local_origin_source_mm=I.tolist(),local_basis_columns=F.tolist(),
        dimensions=dict(head_AF_mm=5.,head_normal_mm=[4.05,8.05],pilot_diameter_mm=1.4,pilot_normal_mm=[3.05,4.05],
            wire_bore_mm=.8,radial_thread='M2x0.4 nominal; 2.1mm display bore is not a tap-drill callout',
            radial_screw_axis_local_z_mm=6.05,nominal_M2_length_mm=3.,flat_point_diameter_mm=.9,
            nominal_tip_local_x_mm=.3,first_touch_only=True),
        parent_204_sources_unchanged=True,internal_wire_and_endpoints_unchanged=True,
        extended_wire_path_source_mm=extended.tolist(),tail_extension_mm=float(np.linalg.norm(tail-path[-1])),
        affected_pairs_checked=count,static_intersections=pairs,static_geometry_pass=not pairs,
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**parent['intended_body_owners'],
            'stop:FJ1507/hex_body':'wire_grip_at_humeral_cap_NOT_verified',
            'stop:FJ1507/M2x3':'stop:FJ1507/hex_body',
            'wire_active:FJ1507/reference':'deformable_active_wire_reference'},
        adopted_as_current=False,stage_pass=False,grip_force_N=None,tightening_torque_Nm=None,actual_drive_range_deg=None,
        not_verified=['body bearing and tool/counterhold access','physical grip, loosening and strand damage',
                      'actual metal threads and manufacture','TPU deformation and driven motion','complete module installation'])
    save(DEST/'assembly_check.json',report)
    print(dict(parts=len(sources),pairs=count,hits=pairs),flush=True)
    if pairs:raise ValueError('Wire stop conflicts with actual source geometry; saved report is not acceptance')


if __name__=='__main__':main()
