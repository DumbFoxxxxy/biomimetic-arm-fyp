"""Independent readback audit of both actual tendon-link packages together.

Keeps free terminal caps distinct from humerus-mounted bodies, records guide
contact probes and cable bend geometry, and forbids using the old rigid-body
sweep as proof of the new flexible mechanism's motion.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_tendon_link as L
import audit_shoulder_retained_joint_fit as A


def minimum_radius(path):
    p=np.array(path);a=p[1:-1]-p[:-2];b=p[2:]-p[1:-1];c=p[2:]-p[:-2]
    cross=np.linalg.norm(np.cross(a,b),axis=1)
    valid=cross>1e-9
    radii=np.linalg.norm(a,axis=1)[valid]*np.linalg.norm(b,axis=1)[valid]*np.linalg.norm(c,axis=1)[valid]/(2*cross[valid])
    return None if not len(radii) else float(min(radii))


def run():
    package_file=L.DEST/'terminal_link_check.json'
    package=json.loads(package_file.read_text(encoding='utf8'))
    if package['script_sha256']!=A.hashed(Path(L.__file__).resolve()):raise ValueError('Stale link generator')
    parent=json.loads(L.SOURCE.read_text(encoding='utf8'))
    if package['parent_assembly_sha256']!=A.hashed(L.SOURCE):raise ValueError('Stale parent assembly')
    sources=dict(parent['sources']);free=[];selection={}
    for head,result in package['results'].items():
        row=result['trials'][result['least_overlap_trial_index']]
        if not result['selected_is_static_clear']:raise ValueError('No clear nominal package: '+head)
        selection[head]=row
        for local,old in row['source_replacements'].items():
            sources[old]=result['sources'][local]
        # Change the cap's identity so no program can silently assume it is
        # still bolted directly to the humerus under its former source key.
        oldcap='cap:'+head+'/distal'
        sources['free_terminal:'+head]=sources.pop(oldcap)
        for local,meta in result['sources'].items():
            if local in row['source_replacements']:continue
            sources['link:'+head+'/'+local]=meta
        free.append(dict(head=head,cap_key='free_terminal:'+head,
            neutral_pin_centre_source_mm=row['cap_pin_centre_source_mm'],
            bone_pin_centre_source_mm=row['bone_pin_centre_source_mm'],
            required='actual flexible muscle deformation/contact and fixed-length routed tendon constraints'))
    solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0]
            for k,r in sources.items()}
    contacts=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:contacts.append(dict(a=a,b=b,volume_mm3=volume,
            category='compliant contact requires review' if a.startswith('tpu:') or b.startswith('tpu:')
            or a.startswith('collar:liner') or b.startswith('collar:liner') else 'rigid/wire fit requires review'))
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes']
    actuation=[]
    for row in routes:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(wire,s)) for k,s in solids.items()]
        actuation.append(dict(head=row['head'],intersections=hits,hits_over_threshold=[h for h in hits if h['volume_mm3']>.001]))
    retention=[];terminal_capture=[]
    for head,row in selection.items():
        r=next(x for x in routes if x['head']==head)
        cap=solids['free_terminal:'+head];muscle=solids['tpu:'+r['module']]
        normal=np.array(r['ports']['distal']['normal_source']);probes=[]
        for sign in [1,-1]:
            first=None
            for step in range(1,21):
                distance=sign*step*.05
                volume=A.I.intersection_volume(cap.translate(normal*distance),muscle)
                if volume>.001:
                    first=dict(translation_mm=distance,intersection_mm3=volume);break
            probes.append(dict(sign=sign,first_interference=first))
        terminal_capture.append(dict(head=head,
            scope='cap normal translations against undeformed TPU pad/hooks; geometric engagement only, NOT clip force, installation or escape proof',
            step_mm=.05,maximum_magnitude_mm=1.,probes=probes))
        collar=json.loads((L.BASE/'interface_check.json').read_text(encoding='utf8'))
        axis=np.array(collar['anchors'][head]['bore_axis_source'])
        bone_axis=axis*row['collar_side'];cap_axis=axis*row['cap_side']
        a=np.array(row['cap_pin_centre_source_mm']);b=np.array(row['bone_pin_centre_source_mm'])
        for end,e in [('muscle',L.eye_data(a,cap_axis,b)),('bone',L.eye_data(b,bone_axis,a))]:
            angles=np.linspace(-np.pi/2,np.pi/2,49)
            arc=e['centre']+2.55*(-np.cos(angles)[:,None]*e['u']+np.sin(angles)[:,None]*e['v'])
            probe=L.wire_tube(arc)
            obstacle=solids['link:'+head+'/'+end+'_groove_ring']
            probes=[]
            for sign in [1,-1]:
                first=None
                for step in range(1,11):
                    distance=sign*step*.05
                    volume=A.I.intersection_volume(probe.translate(e['axis']*distance),obstacle)
                    if volume>.001:
                        first=dict(translation_mm=distance,intersection_mm3=volume);break
                probes.append(dict(sign=sign,first_interference=first))
            retention.append(dict(head=head,end=end,scope='semicircular eye translated along pin axis against actual grooved ring only; NOT grip/load/escape proof',
                                  step_mm=.05,maximum_magnitude_mm=.5,probes=probes))
    # Pin rotations, wire bodies and free caps require different constraints;
    # these labels are declarations of intended assembly, not solved dynamics.
    owners={}
    for key in sources:
        owner='fixed_context'
        if key.startswith('tpu:'):owner='flexible_muscle'
        elif key.startswith('collar:') or key=='bone:humerus' or key.startswith('stop:') or key in ['cap:FJ1468/distal','cap:FJ1507/distal']:
            owner='humerus_ideal_attachment'
        elif key.startswith('free_terminal:'):owner=key
        elif key.startswith('link:'):
            head,part=key[5:].split('/',1)
            owner=('deformable_tendon:'+head if part=='tendon_wire' or 'sleeve' in part else
                   'pin_swivel:'+head+'/'+part if 'groove_ring' in part else
                   'free_terminal:'+head if part.startswith('muscle_') else 'humerus_ideal_attachment')
        for head in selection:
            if key=='stop:'+head+'/distal_stop_blank':owner='free_terminal:'+head
        owners[key]=owner
    report=dict(scope=__doc__,n_instances=len(sources),sources=sources,
        input_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
            [package_file,L.SOURCE,route_file,Path(__file__).resolve()]},
        free_terminal_constraints_required=free,intended_body_owners=owners,
        material_intersections=contacts,old_actuation_wire_checks=actuation,
        finite_guide_probes=retention,finite_terminal_capture_probes=terminal_capture,
        tendon_geometry={h:dict(wire_length_mm=r['nominal_wire_length_mm'],
            centreline_minimum_sampled_radius_mm=minimum_radius(r['wire_path_source_mm']),
            minimum_radius_is_material_approved=False) for h,r in selection.items()},
        prohibited_assumptions=['free terminal rigidly fixed to humerus or world',
            'old straight origin-to-bone-anchor torque proxy is actual routed muscle force',
            'uncrimped sleeves already grip, or nominal M2 envelopes prove fastening/load capacity'],
        not_verified=['new cap equilibrium and coupled deformation, motion and return',
            'crimp grip, wire fatigue, joint/liner/clamp load capacity, mounting and tools',
            'TPU work-axis mapping, nonlinear response and physical drive range'],
        stage_pass=False,actual_drive_range_deg=None)
    (L.DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Together',len(sources),'objects; rigid/wire hits',
        [r for r in contacts if r['volume_mm3']>.001 and r['category'].startswith('rigid')],
        'old actuation wire hits',sum(bool(r['hits_over_threshold']) for r in actuation),flush=True)
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
