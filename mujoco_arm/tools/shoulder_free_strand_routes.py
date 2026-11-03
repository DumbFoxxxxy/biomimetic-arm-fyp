"""Length-matched central strand candidates with real crimp exit tangents.

Retains the nominal 1.6 mm straight exit leads, then uses one symmetric-
handle cubic. This is a declared search family, never a mechanical guide or
proof of a globally shortest feasible route. Installed length is unmeasured.
"""
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.integrate import quad
from scipy.optimize import brentq
import shoulder_free_terminal_constraints as C

def cubic(control,t):
    t=np.asarray(t,float)[...,None];a,b,c,d=control
    return (1-t)**3*a+3*(1-t)**2*t*b+3*(1-t)*t*t*c+t**3*d

def derivative(control,t):
    t=np.asarray(t,float)[...,None];a,b,c,d=control
    return 3*((1-t)**2*(b-a)+2*(1-t)*t*(c-b)+t*t*(d-c))

def arc_length(control):
    value,error=quad(lambda t:float(np.linalg.norm(derivative(control,t))),0.,1.,epsabs=1e-9,epsrel=1e-10,limit=100)
    return float(value),float(error)

def minimum_radius(control):
    a,b,c,d=np.asarray(control,float)
    coefficients=np.array([3*(b-a),6*(c-2*b+a),3*(d-3*c+3*b-a)])
    velocity=[P(coefficients[:,i]) for i in range(3)];accel=[v.deriv() for v in velocity]
    cross=[velocity[1]*accel[2]-velocity[2]*accel[1],velocity[2]*accel[0]-velocity[0]*accel[2],velocity[0]*accel[1]-velocity[1]*accel[0]]
    speed2=sum(v*v for v in velocity);cross2=sum(v*v for v in cross)
    def roots(poly):return [float(z.real) for z in poly.roots() if abs(z.imag)<1e-7 and 0.<z.real<1.]
    min_speed2=min(float(speed2(t)) for t in [0.,1.]+roots(speed2.deriv()))
    if min_speed2<=1e-14:return dict(radius_mm=0.,parameter=None,minimum_speed_squared=min_speed2)
    critical=cross2.deriv()*speed2-3*cross2*speed2.deriv()
    points=[0.,1.]+roots(critical)
    values=[max(0.,float(cross2(t)))/max(float(speed2(t)),1e-30)**3 for t in points]
    i=int(np.argmax(values));radius=None if values[i]<1e-24 else 1./np.sqrt(values[i])
    return dict(radius_mm=None if radius is None else float(radius),parameter=points[i],minimum_speed_squared=min_speed2,
        numerical_critical_points=points,material_bend_limit_verified=False)

def route(model,phiA=0.,phiB=0.,cap_pin=None,cap_rotation=None,humerus_rotation=None):
    Rc=np.eye(3) if cap_rotation is None else cap_rotation;Rh=np.eye(3) if humerus_rotation is None else humerus_rotation
    A,B,_=C.endpoints(model,cap_pin,Rc,Rh,phiA,phiB)
    original=model['active_path'];vA=original[1]-original[0];vB=original[-2]-original[-1]
    leadA=float(np.linalg.norm(vA));leadB=float(np.linalg.norm(vB))
    uA=Rc@C.rotation(model['axisA'],phiA)@(vA/leadA);uB=Rh@C.rotation(model['axisB'],phiB)@(vB/leadB)
    length=model['nominal_active_length'];chord=float(np.linalg.norm(B-A))
    record=dict(phiA_rad=float(phiA),phiB_rad=float(phiB),nominal_length_assumption_mm=length,
        installed_length_measured=False,exit_A_source_mm=A.tolist(),exit_B_source_mm=B.tolist(),
        exit_tangent_A=uA.tolist(),exit_tangent_B_toward_free_segment=uB.tolist(),
        neutral_exit_lead_lengths_mm=[leadA,leadB],exit_chord_mm=chord,
        route_family='source straight leads plus one common-handle clamped cubic')
    if chord>length+1e-6:
        record.update(found=False,reason='exit chord exceeds assumed nominal material length',
            conditional_impossibility_for_prescribed_pose=True)
        return None,record
    p=A+uA*leadA;q=B+uB*leadB;target=length-leadA-leadB;lower=float(np.linalg.norm(q-p))
    record['family_straight_lead_lower_bound_mm']=lower+leadA+leadB
    if lower>=target-1e-9:
        record.update(found=False,reason='retained straight-lead family has no positive-handle length solution',
            conditional_impossibility_for_prescribed_pose=False)
        return None,record
    def controls(h):return np.array([p,p+uA*h,q+uB*h,q])
    upper=max(length,1.)
    while arc_length(controls(upper))[0]<target and upper<length*8:upper*=2
    if arc_length(controls(upper))[0]<target:
        record.update(found=False,reason='finite cubic handle search found no bracket',conditional_impossibility_for_prescribed_pose=False)
        return None,record
    h=brentq(lambda h:arc_length(controls(h))[0]-target,1e-10,upper,xtol=1e-11)
    control=controls(h);arc,error=arc_length(control);curve_radius=minimum_radius(control)
    record.update(handle_mm=float(h),control_points_source_mm=control.tolist(),cubic_arc_length_mm=arc,
        arc_quadrature_error_estimate_mm=error,arc_length_error_mm=arc+leadA+leadB-length,
        minimum_cubic_radius=curve_radius)
    # A radius below the wire radius makes this swept-tube construction
    # locally singular; this is a geometry guard, not a steel bend criterion.
    if curve_radius['radius_mm'] is not None and curve_radius['radius_mm']<=.300001:
        record.update(found=False,reason='candidate curve radius not above 0.3mm wire radius',conditional_impossibility_for_prescribed_pose=False)
        return None,record
    n=33
    while True:
        path=C.deduplicate(np.vstack([A,cubic(control,np.linspace(0,1,n)),B]))
        deficit=length-C.polyline_length(path)
        if deficit<=1e-4 or n>=4097:break
        n=2*n-1
    if deficit>1e-4:
        record.update(found=False,reason='curve sampling length tolerance not met',conditional_impossibility_for_prescribed_pose=False)
        return None,record
    record.update(found=True,reason=None,conditional_impossibility_for_prescribed_pose=False,
        sampled_path_source_mm=path.tolist(),polyline_length_deficit_mm=deficit,sample_count=len(path),
        minimum_radius_is_NOT_a_manufacturer_approved_bend_radius=True)
    return path,record
