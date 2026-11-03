"""Source-based shoulder tendon endpoints and tension-only chord constraint.

The nominal free-strand material length is extracted BETWEEN the two crimp
outer exits, not from the full folded wire (which includes eyes and tails).
No-slip crimps, inextensibility, and a straight unobstructed taut strand are
conditional model assumptions, not verified physical properties.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import design_shoulder_tendon_link as L

def rotation(axis,angle_rad):
    a=np.asarray(axis,float);a=a/np.linalg.norm(a)
    K=np.array([[0.,-a[2],a[1]],[a[2],0.,-a[0]],[-a[1],a[0],0.]])
    return np.eye(3)+np.sin(angle_rad)*K+(1-np.cos(angle_rad))*(K@K)

def polyline_length(path):return float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum())

def locate(path,point):
    a=path[:-1];v=path[1:]-a;den=np.sum(v*v,axis=1)
    t=np.clip(np.sum((point-a)*v,axis=1)/np.maximum(den,1e-30),0.,1.)
    error=np.linalg.norm(a+t[:,None]*v-point,axis=1);i=int(np.argmin(error))
    if error[i]>1e-6:raise ValueError('Crimp exit not on recorded centreline')
    return i,float(t[i]),float(error[i])

def deduplicate(p):
    p=np.asarray(p,float)
    return p[np.r_[True,np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-9]]

def source_model():
    selection_file=A.ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=json.loads(selection_file.read_text(encoding='utf8'))
    assembly_file=A.ROOT/selection['assembly_report'];assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    if A.hashed(assembly_file)!=selection['assembly_sha256']:raise ValueError('Stale selection')
    package_file=L.DEST/'terminal_link_check.json';package=json.loads(package_file.read_text(encoding='utf8'))
    if package['script_sha256']!=A.hashed(Path(L.__file__).resolve()):raise ValueError('Changed link generator')
    collar_file=L.BASE/'interface_check.json';collar=json.loads(collar_file.read_text(encoding='utf8'))
    heads=[]
    for required in assembly['free_terminal_constraints_required']:
        head=required['head'];entry=package['results'][head];trial=entry['trials'][entry['least_overlap_trial_index']]
        path=np.array(trial['wire_path_source_mm']);pinA=np.array(trial['cap_pin_centre_source_mm']);pinB=np.array(trial['bone_pin_centre_source_mm'])
        if max(np.linalg.norm(pinA-required['neutral_pin_centre_source_mm']),np.linalg.norm(pinB-required['bone_pin_centre_source_mm']))>1e-8:raise ValueError('Changed pin locations')
        axis=np.array(collar['anchors'][head]['bore_axis_source']);axisA=axis*trial['cap_side'];axisB=axis*trial['collar_side']
        eyeA=L.eye_data(pinA,axisA,pinB);eyeB=L.eye_data(pinB,axisB,pinA)
        rebuilt=L.folded_wire(eyeA,eyeB)
        if rebuilt.shape!=path.shape or np.max(abs(rebuilt-path))>1e-8:raise ValueError('Recorded routing does not reproduce')
        exitA=eyeA['centre']+eyeA['u']*8.9+eyeA['v']*.45
        exitB=eyeB['centre']+eyeB['u']*8.9+eyeB['v']*.45
        ia,ta,ea=locate(path,exitA);ib,tb,eb=locate(path,exitB)
        if ia>=ib:raise ValueError('Unexpected wire ordering')
        active=deduplicate(np.vstack([exitA,path[ia+1:ib+1],exitB]))
        prefix=deduplicate(np.vstack([path[:ia+1],exitA]));suffix=deduplicate(np.vstack([exitB,path[ib+1:]]))
        verified_keys=[]
        for local,meta in entry['sources'].items():
            if local in ['muscle_TPU_with_hooks']:continue
            key=('free_terminal:'+head if local=='free_terminal_cap' else
                 trial['source_replacements'][local] if local in trial['source_replacements'] else 'link:'+head+'/'+local)
            if key not in assembly['sources'] or assembly['sources'][key]['sha256']!=meta['sha256']:raise ValueError('Changed tendon hardware '+key)
            if A.hashed(A.ROOT/meta['file'])!=meta['sha256']:raise ValueError('Changed hardware file')
            verified_keys.append(key)
        heads.append(dict(head=head,pinA=pinA,pinB=pinB,axisA=axisA,axisB=axisB,exitA=exitA,exitB=exitB,
            full_path=path,active_path=active,prefix=prefix,suffix=suffix,
            nominal_active_length=polyline_length(active),nominal_full_length=polyline_length(path),
            exit_location_errors_mm=[ea,eb],exit_polyline_parameters=[[ia,ta],[ib,tb]],verified_hardware_keys=verified_keys))
    return assembly_file,assembly,heads,dict(selection_sha256=A.hashed(selection_file),package_sha256=A.hashed(package_file),collar_sha256=A.hashed(collar_file))

def endpoints(model,cap_pin=None,cap_rotation=None,humerus_rotation=None,phiA=0.,phiB=0.):
    x=model['pinA'] if cap_pin is None else np.asarray(cap_pin,float)
    Rc=np.eye(3) if cap_rotation is None else cap_rotation;Rh=np.eye(3) if humerus_rotation is None else humerus_rotation
    G=np.asarray(A.frames.joint_centres('source')['GH'])
    offsetA=rotation(model['axisA'],phiA)@(model['exitA']-model['pinA'])
    offsetB=rotation(model['axisB'],phiB)@(model['exitB']-model['pinB'])
    a=x+Rc@offsetA;bpin=G+Rh@(model['pinB']-G);b=bpin+Rh@offsetB
    return a,b,bpin

def chord_and_jacobian(model,cap_pin=None,cap_rotation=None,humerus_rotation=None,phiA=0.,phiB=0.):
    Rc=np.eye(3) if cap_rotation is None else cap_rotation;Rh=np.eye(3) if humerus_rotation is None else humerus_rotation
    x=model['pinA'] if cap_pin is None else np.asarray(cap_pin,float)
    a,b,bpin=endpoints(model,x,Rc,Rh,phiA,phiB)
    length=float(np.linalg.norm(a-b))
    if length<1e-9:raise ValueError('Degenerate chord')
    e=(a-b)/length;ra=a-x;rb=b-bpin
    J=np.r_[e,np.cross(ra,e),e@np.cross(Rc@model['axisA'],ra),-e@np.cross(Rh@model['axisB'],rb)]
    return length,J

def evaluate_installed_strand(model,installed_free_strand_length_mm,*,length_provenance,
        clear_routed_length_mm=None,straight_chord_clear_verified=False,cap_pin=None,
        cap_rotation=None,humerus_rotation=None,phiA=0.,phiB=0.):
    """Evaluate a conditional length gap, never solve force or a cap pose.

    Caller supplies a measured or explicitly assumed INSTALLED free length.
    A shortest feasible routed length is also required; a chord is used only
    after explicit clearance verification for that exact prescribed pose.
    """
    if installed_free_strand_length_mm is None:
        raise ValueError('Installed free-strand length is unknown; full folded-wire length is not a substitute')
    installed=float(installed_free_strand_length_mm)
    if not np.isfinite(installed) or installed<=0:raise ValueError('Need positive finite installed free-strand length')
    if length_provenance not in ['measured','explicit_assumption']:
        raise ValueError('Declare measured or explicit_assumption length provenance')
    chord,_=chord_and_jacobian(model,cap_pin,cap_rotation,humerus_rotation,phiA,phiB)
    if clear_routed_length_mm is None:
        if not straight_chord_clear_verified:raise ValueError('Need an obstacle-cleared route length or verified clear chord at this pose')
        routed=chord;route_kind='explicitly clearance-verified straight chord'
    else:
        routed=float(clear_routed_length_mm);route_kind='caller-provided shortest feasible routed length'
        if not np.isfinite(routed) or routed<chord-1e-6:raise ValueError('Feasible routed length cannot be shorter than the exit chord')
    gap=installed-routed
    branch=('positive geometric slack under conditional ideal assumptions' if gap>1e-6 else
            'prescribed state incompatible with inextensible no-slip assumption' if gap < -1e-6 else
            'taut or zero tension; force requires equilibrium')
    return dict(installed_free_strand_length_mm=installed,length_provenance=length_provenance,
        neutral_design_length_is_not_automatically_installed=True,chord_length_mm=chord,
        routed_length_mm=routed,route_kind=route_kind,length_gap_mm=gap,conditional_branch=branch,
        tension_N=None,cap_pose_solved=False,actual_drive_range_deg=None,stage_pass=False)
