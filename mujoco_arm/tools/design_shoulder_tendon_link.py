"""Actual nominal free-terminal/tendon hardware for the shaft-anchor candidate.

The folded wire has two eyes and two UNCRIMPED double-barrel sleeves. Groove
rings and M2 pin envelopes provide actual geometric guidance; no unsupported
world-fixed guide or rigidly humerus-mounted free cap is imposed. This is a
neutral-assembly design, not measured grip, flexible equilibrium or motion.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_shaft_collar import rod
from build_shoulder_isolated_retained_joint import export_checked
from build_shoulder_unified_muscle import local_export

BASE=A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60'
SOURCE=BASE/'deltoid_relief/assembly_check.json'
DEST=BASE/'tendon_links'


def basis(axis):
    z=A.D.unit(axis);helper=np.eye(3)[np.argmin(abs(z))]
    x=A.D.unit(np.cross(helper,z));y=np.cross(z,x)
    return np.stack([x,y,z],axis=1)


def mapped(solid,axis,centre):
    return solid.transform(np.c_[basis(axis),np.array(centre)])


def grooved_ring(centre,axis):
    """R2.55 wire path, 0.35-radius groove, actual M2-bearing bore."""
    mf=A.I._manifold()
    ring=mf.Manifold.cylinder(1.,2.75,2.75,64).translate([0,0,-.5])
    ring=ring-mf.Manifold.cylinder(2.,1.05,1.05,48).translate([0,0,-1.])
    angle=np.linspace(0,2*np.pi,33)[:-1]
    profile=np.c_[2.55+.35*np.cos(angle),.35*np.sin(angle)]
    groove=mf.Manifold.revolve(mf.CrossSection([profile]),64)
    return mapped(ring-groove,axis,centre)


def m2_pin(centre,axis):
    """M2x10 nominal screw/nut/washer envelopes, no thread or locking proof."""
    mf=A.I._manifold();axis=A.D.unit(axis);F=basis(axis)
    profile=np.array([[0,-6.6],[1.,-6.6],[1.,3.4],[1.9,3.4],[1.9,5.4],[0,5.4]])
    bolt=mapped(mf.Manifold.revolve(mf.CrossSection([profile]),64),axis,centre)
    nut=mf.Manifold.cylinder(1.6,4/np.sqrt(3),4/np.sqrt(3),6).translate([0,0,-3.8])
    nut=nut-mf.Manifold.cylinder(3.,1.1,1.1,48).translate([0,0,-4.5])
    parts={'bolt':bolt,'nut':mapped(nut,axis,centre)}
    for name,start,end in [('inner_washer',-2.1,-1.6),('outer_washer',2.8,3.3)]:
        parts[name]=(rod(centre+axis*start,centre+axis*end,2.)-
                     rod(centre+axis*(start-.2),centre+axis*(end+.2),1.1))
    return parts


def bezier(p0,p1,p2,p3,n=12):
    t=np.linspace(0,1,n)[:,None]
    return (1-t)**3*p0+3*(1-t)**2*t*p1+3*(1-t)*t*t*p2+t**3*p3


def eye_data(pin_centre,axis,toward):
    axis=A.D.unit(axis);centre=np.array(pin_centre)+axis*2.2
    delta=np.array(toward)-centre;u=A.D.unit(delta-axis*(delta@axis));v=np.cross(axis,u)
    R=2.55;offset=.45
    return dict(centre=centre,u=u,v=v,axis=axis,
        minus=centre-v*R,plus=centre+v*R,
        tail=centre+u*10-v*offset,qtail=centre+u*6-v*offset,
        qmain=centre+u*6+v*offset)


def folded_wire(e0,e1):
    def arc(e):
        a=np.linspace(-np.pi/2,np.pi/2,49)
        return e['centre']+2.55*(-np.cos(a)[:,None]*e['u']+np.sin(a)[:,None]*e['v'])
    pieces=[np.array([e0['tail'],e0['qtail']]),
        bezier(e0['qtail'],e0['qtail']-e0['u']*2,e0['minus']+e0['u']*2,e0['minus']),
        arc(e0),bezier(e0['plus'],e0['plus']+e0['u']*2,e0['qmain']-e0['u']*2,e0['qmain'])]
    lead0=e0['qmain']+e0['u']*4.5;lead1=e1['qmain']+e1['u']*4.5
    length=np.linalg.norm(lead1-lead0)
    pieces.extend([np.array([e0['qmain'],lead0]),
        bezier(lead0,lead0+e0['u']*length/3,lead1+e1['u']*length/3,lead1,32),
        np.array([lead1,e1['qmain']])])
    pieces.extend([bezier(e1['qmain'],e1['qmain']-e1['u']*2,e1['plus']+e1['u']*2,e1['plus']),
        arc(e1)[::-1],bezier(e1['minus'],e1['minus']+e1['u']*2,e1['qtail']-e1['u']*2,e1['qtail']),
        np.array([e1['qtail'],e1['tail']])])
    path=np.vstack(pieces);keep=np.r_[True,np.linalg.norm(np.diff(path,axis=0),axis=1)>1e-8]
    return path[keep]


def sleeve(e):
    parts=[];holes=[]
    for sign in [-1,1]:
        offset=e['v']*.45*sign
        parts.append(rod(e['centre']+e['u']*6.4+offset,e['centre']+e['u']*8.9+offset,1.))
        holes.append(rod(e['centre']+e['u']*6+offset,e['centre']+e['u']*9.3+offset,.4))
    return (parts[0]+parts[1])-holes[0]-holes[1]


def wire_tube(path,radius=.3,sides=16):
    """One indexed tube along the continuous cable, without CSG seam debris."""
    tangents=np.vstack([path[1]-path[0],path[2:]-path[:-2],path[-1]-path[-2]])
    tangents/=np.linalg.norm(tangents,axis=1)[:,None]
    theta=np.arange(sides)*2*np.pi/sides;vertices=[]
    x=basis(tangents[0])[:,0]
    for point,t in zip(path,tangents):
        x=A.D.unit(x-t*(x@t));y=np.cross(t,x)
        vertices.extend(point+radius*(np.cos(theta)[:,None]*x+np.sin(theta)[:,None]*y))
    vertices.extend([path[0],path[-1]]);faces=[];count=len(path)
    for i in range(count-1):
        for j in range(sides):
            a=i*sides+j;b=i*sides+(j+1)%sides;c=(i+1)*sides+j;d=(i+1)*sides+(j+1)%sides
            faces.extend([[a,b,c],[c,b,d]])
    for j in range(sides):
        faces.append([count*sides,(j+1)%sides,j])
        faces.append([count*sides+1,(count-1)*sides+j,(count-1)*sides+(j+1)%sides])
    return A.I._manifold().Manifold(A.I._manifold().Mesh64(np.array(vertices),np.array(faces,dtype=np.uint32)))


def clearance_union(solids,clearance=.1):
    pieces=[]
    for solid in solids:
        pieces.append(solid)
        pieces.extend(solid.translate(np.eye(3)[axis]*sign*clearance)
                      for axis in range(3) for sign in [-1,1])
    return A.I._manifold().Manifold.batch_boolean(pieces,A.I._manifold().OpType.Add)


def candidate(head,side_cap,side_collar,assembly,collar,routes,source_solids):
    row=next(r for r in routes if r['head']==head)
    port=row['ports']['distal'];spec=row['distal_terminal_override']
    n=np.array(port['normal_source']);inner=np.array(port['inner_face_source_mm'])
    bolt=np.array(spec['bolt_inner_source_mm']);axis=np.array(collar['anchors'][head]['bore_axis_source'])
    cap_axis=axis*side_cap;bone_axis=axis*side_collar
    cap_pin=bolt+n*6.5;bone_pin=np.array(collar['anchors'][head]['eye_centre_source_mm'])
    oldcap=source_solids['cap:'+head+'/distal']
    # Fill the unused M3 mounting lobe and make an actual integral tendon ear.
    ear=rod(cap_pin-cap_axis*1.5,cap_pin+cap_axis*1.5,3.2)
    cap=oldcap+rod(bolt+n*.15,cap_pin-n*1.2,2.3)+ear
    cap=cap-rod(cap_pin-cap_axis*5,cap_pin+cap_axis*5,1.1)
    e0=eye_data(cap_pin,cap_axis,bone_pin);e1=eye_data(bone_pin,bone_axis,cap_pin)
    path=folded_wire(e0,e1)
    parts={'free_terminal_cap':cap,'tendon_wire':wire_tube(path)}
    for end,pin,paxis,e in [('muscle',cap_pin,cap_axis,e0),('bone',bone_pin,bone_axis,e1)]:
        parts.update({end+'_'+k:v for k,v in m2_pin(pin,paxis).items()})
        parts[end+'_groove_ring']=grooved_ring(e['centre'],paxis)
        parts[end+'_uncrimped_sleeve']=sleeve(e)
    jaw_key='collar:jaw_'+('upper' if head=='FJ1500' else 'lower')
    oldjaw=source_solids[jaw_key]
    cap_cut=clearance_union([v for k,v in parts.items()
        if k in ['muscle_nut','muscle_inner_washer','muscle_outer_washer','muscle_groove_ring','tendon_wire']])
    jaw_cut=clearance_union([v for k,v in parts.items()
        if k in ['bone_nut','bone_inner_washer','bone_outer_washer','bone_groove_ring','tendon_wire']])
    parts['free_terminal_cap']=cap-cap_cut
    parts['collar_jaw_revised']=oldjaw-jaw_cut
    # The old bone-side M3 lobe is now an integral tendon ear. Replace its
    # lost TPU retention with actual flexible hooks rather than relying only
    # on an unsupported claim that a loose cap remains attached to the pad.
    tangent=np.array(spec['tangent_source']);tangent=A.D.unit(tangent-n*(tangent@n))
    lateral=np.cross(n,tangent);F=np.stack([tangent,lateral,n],axis=1)
    def local_box(lo,hi):
        return A.I._manifold().Manifold.cube(np.array(hi)-lo).translate(lo).transform(np.c_[F,inner])
    hooks=[]
    for sign in [1,-1]:
        yl,yh=(3.4,5.) if sign>0 else (-5.,-3.4)
        lip_l,lip_h=(2.6,5.) if sign>0 else (-5.,-2.6)
        hooks.extend([local_box([-1.,yl,-2.],[1.,yh,4.]),
                      local_box([-1.,lip_l,3.2],[1.,lip_h,4.])])
    muscle_key='tpu:'+row['module'];old_muscle=source_solids[muscle_key]
    parts['muscle_TPU_with_hooks']=old_muscle+A.I._manifold().Manifold.batch_boolean(hooks,A.I._manifold().OpType.Add)
    parts['muscle_TPU_with_hooks']=parts['muscle_TPU_with_hooks']-A.D.outer_envelope(
        np.array(row['full_wire_path_source_mm']),.86)
    bearing_protection={}
    for name,point,paxis,old,new in [('muscle',cap_pin,cap_axis,cap,parts['free_terminal_cap']),
                                   ('bone',bone_pin,bone_axis,oldjaw,parts['collar_jaw_revised'])]:
        bearing=(rod(point-paxis*1.5,point+paxis*1.5,3.2)-
                 rod(point-paxis*2,point+paxis*2,1.1))
        bearing_protection[name]=A.I.intersection_volume(old-new,bearing)
    replacement='cap:'+head+'/distal'
    conflicts=[]
    for a,b in itertools.combinations(parts,2):
        volume=A.I.intersection_volume(parts[a],parts[b])
        if volume>.001:conflicts.append(dict(a=a,b=b,volume_mm3=volume))
    external=[]
    for a,solid in parts.items():
        for b,other in source_solids.items():
            if b in [replacement,jaw_key,muscle_key]:continue
            volume=A.I.intersection_volume(solid,other)
            if volume>.001:
                added_overlap=(A.I.intersection_volume(solid-old_muscle,other)
                               if a=='muscle_TPU_with_hooks' else volume)
                external.append(dict(part=a,reference=b,volume_mm3=volume,
                    new_material_overlap_mm3=added_overlap,
                    category='TPU contact review' if a=='muscle_TPU_with_hooks' or b.startswith('tpu:') or b.startswith('collar:liner') else 'rigid/wire conflict'))
    # Internal actuation wire must retain its old port/stop clearance.
    oldwire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
    actuation_hits=[dict(part=k,volume_mm3=A.I.intersection_volume(oldwire,s)) for k,s in parts.items()]
    failed_internal=[r for r in actuation_hits if r['volume_mm3']>.001]
    score=sum(r['volume_mm3'] for r in conflicts)+sum(r['new_material_overlap_mm3'] for r in external)+sum(r['volume_mm3'] for r in failed_internal)
    return parts,dict(head=head,cap_side=side_cap,collar_side=side_collar,
        cap_pin_centre_source_mm=cap_pin.tolist(),bone_pin_centre_source_mm=bone_pin.tolist(),
        free_cap_is_not_humerus_rigidly_attached=True,wire_path_source_mm=path.tolist(),
        nominal_wire_length_mm=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()),
        internal_part_conflicts=conflicts,external_conflicts=external,
        original_actuation_wire_intersections=actuation_hits,
        all_static_material_checks_clear=(not conflicts and not any(x['new_material_overlap_mm3']>.001 for x in external) and not failed_internal
            and all(v<=1e-6 for v in bearing_protection.values())),
        relief_removed_mm3=dict(cap=float((cap-parts['free_terminal_cap']).volume()),
                                collar_jaw=float((oldjaw-parts['collar_jaw_revised']).volume())),
        protected_bearing_material_removed_mm3=bearing_protection,
        TPU_hooks=dict(added_material_mm3=float((parts['muscle_TPU_with_hooks']-old_muscle).volume()),
            removed_from_old_TPU_mm3=float((old_muscle-parts['muscle_TPU_with_hooks']).volume()),
            preserved_actuation_lumen_radius_mm=.86,
            lip_clearance_mm=.2,lip_overlap_mm=.6,clip_stem_cross_section_mm=[2.,1.6],
            snap_deflection_force_and_fatigue_verified=False),
        source_replacements={'free_terminal_cap':replacement,'collar_jaw_revised':jaw_key,'muscle_TPU_with_hooks':muscle_key},
        overlap_score_mm3=score,stage_pass=False)


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    assembly=json.loads(SOURCE.read_text(encoding='utf8'))
    collar_file=BASE/'interface_check.json';collar=json.loads(collar_file.read_text(encoding='utf8'))
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes']
    solids={k:A.load_checked(A.ROOT/v['file'],v['translation_source_mm'],v['sha256'])[0]
            for k,v in assembly['sources'].items()}
    result={}
    for head in ['FJ1500','FJ1504']:
        trials=[];built=[]
        for a,b in itertools.product([1,-1],repeat=2):
            parts,trial=candidate(head,a,b,assembly,collar,routes,solids)
            trials.append(trial);built.append(parts)
            print(head,a,b,'overlap',round(trial['overlap_score_mm3'],5),
                'internal',[(x['a'],x['b'],round(x['volume_mm3'],4)) for x in trial['internal_part_conflicts']],
                'external',[(x['part'],x['reference'],round(x['volume_mm3'],4)) for x in trial['external_conflicts']],flush=True)
        index=min(range(len(trials)),key=lambda i:trials[i]['overlap_score_mm3'])
        selected=trials[index];folder=DEST/head;folder.mkdir(exist_ok=True)
        exported={}
        for key,solid in built[index].items():
            if key=='muscle_TPU_with_hooks':
                original=assembly['sources'][selected['source_replacements'][key]]
                origin=np.array(original['translation_source_mm'] or A.frames.SHIFT)
                actual,exported[key]=local_export(solid,folder/(key+'.stl'),origin)
                if not exported[key]['mesh_checks_pass']:raise ValueError('TPU hook mesh fails readback: '+head)
            else:
                actual,exported[key]=export_checked(solid.as_original().simplify(.001),folder/(key+'.stl'))
            exported[key]['kind']='TPU_lattice_with_terminal_hooks' if key=='muscle_TPU_with_hooks' else 'nominal_wire_NOT_print' if key=='tendon_wire' else 'uncrimped_sleeve' if 'sleeve' in key else 'rigid_candidate'
        result[head]=dict(trials=trials,least_overlap_trial_index=index,
            selected_is_static_clear=selected['all_static_material_checks_clear'],sources=exported)
    report=dict(scope=__doc__,parent_assembly_sha256=A.hashed(SOURCE),collar_sha256=A.hashed(collar_file),
        route_sha256=A.hashed(route_file),script_sha256=A.hashed(Path(__file__).resolve()),results=result,
        not_verified=['pressing/crimp grip and fatigue; sleeves are uncrimped blanks',
            'M2 thread/locking, printed clearances, wear and load capacity',
            'flexible muscle/cap equilibrium and full motion with new tendons',
            'new free terminals may not be treated as humerus-mounted rigid caps'],
        stage_pass=False,actual_drive_range_deg=None)
    (DEST/'terminal_link_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
