"""Single-M3 compact wire stop plate, geometry candidate only.

Nominal uncrimped stop sleeve is captured against a 0.4mm shoulder. No
claim about grip, preload, anti-rotation, bone fixation or strength.
"""
import numpy as np
import design_shoulder_port_paths as D
import build_shoulder_tendon_segments as T
import build_shoulder_ends as BE
import interference as I


def housing(inner,normal,tangent,radius=3.2,thickness=3.):
    inner=np.asarray(inner,float);n=D.unit(normal)
    tangent=np.asarray(tangent,float);t=D.unit(tangent-n*(tangent@n));b=np.cross(n,t)
    pitch=5.;bolt=inner+t*pitch
    cap=BE.disc(inner,n,thickness,radius=radius,holes=False)+BE.disc(bolt,n,thickness,radius=radius,holes=False)
    bridge=I._manifold().Manifold.cube([pitch,2*radius,thickness]).translate([0,-radius,0])
    cap=cap+bridge.transform(np.c_[np.stack([t,b,n],axis=1),inner])
    return cap


def build(inner,normal,tangent):
    inner=np.asarray(inner,float);n=D.unit(normal)
    tangent=np.asarray(tangent,float);t=D.unit(tangent-n*(tangent@n))
    thickness=3.;radius=3.2;pitch=5.;outer=inner+n*thickness;bolt=inner+t*pitch
    cap=housing(inner,n,t)
    cap=cap-T.rod(inner-n*3,outer+n*3,.86)-T.rod(bolt-n*3,bolt+n*6,1.7)
    cap=cap-T.rod(outer-n*2.6,outer+n,1.4)
    stop=T.rod(outer-n*2.55,outer-n*.05,1.25)-T.rod(outer-n*3,outer+n*.2,.45)
    return cap,stop,dict(inner_face_source_mm=inner.tolist(),outer_face_source_mm=outer.tolist(),
        normal_source=n.tolist(),tangent_source=t.tolist(),bolt_inner_source_mm=bolt.tolist(),
        thickness_mm=thickness,outer_radius_mm=radius,bolt_pitch_from_wire_mm=pitch,
        bolt_hole_radius_mm=1.7,wire_bore_radius_mm=.86,stop_counterbore_radius_mm=1.4,
        stop_counterbore_depth_mm=2.6,stop_outer_radius_mm=1.25,stop_length_mm=2.5,
        nominal_shoulder_thickness_mm=.4,nominal_axial_stop_gap_mm=.05,
        fixation_verified=False,crimp_grip_verified=False)
