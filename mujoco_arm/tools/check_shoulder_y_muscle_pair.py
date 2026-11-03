"""Independent current-assembly, paper-frame and interface probes of Y bellies.

Finite normal contact and static wire checks do not validate deformation,
material force/strain, printing or the free-terminal/tendon equilibrium.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_y_muscle_pair import DEST
from build_shoulder_shaft_collar import rod
from check_shoulder_front_liner_keepers import is_compliant


def run():
    source=DEST/'candidate_generation.json';generation=json.loads(source.read_text(encoding='utf8'))
    if not generation['both_muscles_have_geometry_candidate']:raise ValueError('Need both geometry candidates')
    if generation['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_y_muscle_pair.py'):
        raise ValueError('Changed Y muscle generator')
    if not all(A.hashed(A.ROOT/p)==v for p,v in generation['generation_dependency_sha256'].items()):
        raise ValueError('Changed generation dependency')
    parent_file=A.ROOT/generation['parent_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=generation['parent_sha256']:raise ValueError('Changed baseline')
    S={};census={}
    for key,row in generation['sources'].items():
        S[key],_,meta=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
        positive=sum(x.volume()>1e-6 for x in S[key].decompose())
        if positive!=1:raise ValueError('Disconnected material '+key)
        census[key]=dict(boundary=meta['boundary_edges'],nonmanifold=meta['nonmanifold_edges'],positive_components=positive)
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=json.loads(routefile.read_text(encoding='utf8'))['routes']
    overlaps=[]
    for a,b in itertools.combinations(S,2):
        v=A.I.intersection_volume(S[a],S[b])
        if v>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=v,compliant=is_compliant(a) or is_compliant(b)))
    wires=[]
    for r in routes:
        wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)
        hits={k:A.I.intersection_volume(wire,s) for k,s in S.items()}
        wires.append(dict(head=r['head'],volumes_mm3=hits,hits_over_threshold=[dict(part=k,volume_mm3=v) for k,v in hits.items() if v>.001]))
    fixation_file=A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/origin_fixation_review/fixation_audit.json'
    fixation=json.loads(fixation_file.read_text(encoding='utf8'));muscles=[]
    for row in generation['attempts']:
        if not row['selected']:continue
        head=row['head'];key='tpu:'+row['module'];r=next(r for r in routes if r['head']==head)
        oldrow=parent['sources'][key];old=A.load_checked(A.ROOT/oldrow['file'],oldrow['translation_source_mm'],oldrow['sha256'])[0]
        P=np.array(r['core_path_source_mm']);dist=np.linalg.norm(np.diff(P,axis=0),axis=1);s=np.r_[0,np.cumsum(dist)]
        station_checks=[];frames=row['transported_frames']['stations']
        for k,frame in enumerate(frames):
            target=s[-1]*k/(len(frames)-1);i=min(max(0,int(np.searchsorted(s,target,side='left'))-1),len(P)-2)
            tangent=(P[i+1]-P[i])/dist[i];expected=P[i]+tangent*(target-s[i]);F=np.array(frame['paper_basis_columns_source'])
            station_checks.append(dict(station=k,centre_error_mm=float(np.linalg.norm(np.array(frame['centre_source_mm'])-expected)),
                                       orthogonality_error=float(np.max(abs(F.T@F-np.eye(3)))),determinant=float(np.linalg.det(F)),
                                       contraction_Y_dot_path_tangent=float(F[:,1]@tangent)))
        graph=row['transported_frames'];node_map={tuple(v['index']):np.array(v['source_mm']) for v in graph['nodes']}
        slab_checks={k:dict(index=k,cell_count=0,cells_with_three_supported_branches=0) for k in range(len(frames)-1)}
        cell_checks=[]
        for cell in graph['cells']:
            centre=np.array(cell['centre_source_mm']);supported=0;branch_fractions=[]
            for corner in cell['corner_indices']:
                q=node_map[tuple(corner)];fractions=[]
                for fraction in [.35,.65]:
                    point=centre+(q-centre)*fraction
                    probe=A.I._manifold().Manifold.sphere(.15,16).translate(point)
                    fractions.append(A.I.intersection_volume(probe,S[key])/probe.volume())
                branch_fractions.append(max(fractions))
                if max(fractions)>.8:supported+=1
            slab_checks[cell['index'][2]]['cell_count']+=1
            if supported>=3:slab_checks[cell['index'][2]]['cells_with_three_supported_branches']+=1
            cell_checks.append(dict(index=cell['index'],supported_branch_count=supported,
                                    maximum_material_fraction_per_branch=branch_fractions))
        pitches={'X':[],'Z':[]};nx,nz=graph['cell_count_X_Z']
        surface_checks=[];surface_slabs={i:0 for i in range(len(frames)-1)}
        for face in graph['surface_faces']:
            centre=np.array(face['centre_source_mm']);branches=0
            for corner in face['corner_indices']:
                q=node_map[tuple(corner)];ratios=[]
                for fraction in [.4,.7]:
                    probe=A.I._manifold().Manifold.sphere(.1,16).translate(centre+(q-centre)*fraction)
                    ratios.append(A.I.intersection_volume(probe,S[key])/probe.volume())
                if max(ratios)>.8:branches+=1
            if branches>=2:surface_slabs[face['axial_slab']]+=1
            surface_checks.append(dict(axial_slab=face['axial_slab'],perimeter_edge=face['perimeter_edge'],supported_branches=branches))
        for k,f in enumerate(frames):
            F=np.array(f['paper_basis_columns_source'])
            for i in range(nx+1):
                for j in range(nz+1):
                    p=node_map[i,j,k]
                    if i<nx:pitches['X'].append(abs(float((node_map[i+1,j,k]-p)@F[:,0])))
                    if j<nz:pitches['Z'].append(abs(float((node_map[i,j+1,k]-p)@F[:,2])))
        n=np.array(r['ports']['origin']['normal_source']);cap='cap:'+head+'/origin';front='origin:'+head+'/front_support'
        contacts=[]
        for part,direction in [(cap,-n),(front,n)]:
            contacts.append(dict(part=part,translation_source_mm=(direction*.05).tolist(),
                                 old_gap_capped_1mm=A.I.gap(old,S[part],1.),new_gap_capped_1mm=A.I.gap(S[key],S[part],1.),
                                 old_0_05mm_probe_overlap_mm3=A.I.intersection_volume(S[part].translate(direction*.05),old),
                                 new_0_05mm_probe_overlap_mm3=A.I.intersection_volume(S[part].translate(direction*.05),S[key])))
        holes=[];pad_height=2.5 if head=='FJ1500' else 4.
        holespec=next(v for v in fixation['rows'] if v['head']==head)
        for hole in holespec['holes']:
            q=np.array(hole['axis_point_source_mm']);normal=np.array(hole['normal_source'])
            probe=rod(q-normal*(pad_height+1.),q+normal*.5,1.5)
            holes.append(dict(sign=hole['sign'],M3_shaft_local_pad_overlap_mm3=A.I.intersection_volume(probe,S[key])))
        cap=S['free_terminal:'+head];normal=np.array(r['ports']['distal']['normal_source']);capture=[]
        for sign in [-1,1]:
            first=None
            for distance in np.arange(.05,1.001,.05):
                value=A.I.intersection_volume(cap.translate(normal*(distance*sign)),S[key])
                if value>.001:first=dict(distance_mm=round(float(distance*sign),6),volume_mm3=value);break
            capture.append(dict(sign=sign,first_interference=first))
        full=np.array(r['full_wire_path_source_mm'])
        lumen=A.D.outer_envelope(full,.86);guard=A.D.outer_envelope(P,1.61)-lumen
        ends=[]
        for port in r['ports'].values():
            point=np.array(port['inner_face_source_mm']);normal=np.array(port['normal_source'])
            ends.append(rod(point-normal*8.,point+normal*3.,12.))
        exclusion=ends[0]+ends[1];guard-=exclusion
        old_missing=guard-old;new_missing=guard-S[key]
        stabilized_guard=(A.D.outer_envelope(P,1.612)-A.D.outer_envelope(full,.862))-exclusion
        stable_missing=stabilized_guard-S[key]
        muscles.append(dict(module=row['module'],head=head,nominal_cell_xyz_mm=row['cell_xyz_mm'],
                            selected_cells=row['cells'],actual_station_spacing_mm=row['transported_frames']['station_spacing_mm'],
                            transverse_pitch_ranges_mm={axis:[min(v),max(v)] for axis,v in pitches.items()},
                            axial_slab_coverage=list(slab_checks.values()),cell_branch_material_checks=cell_checks,
                            all_axial_slabs_have_supported_cells=all(v['cells_with_three_supported_branches']>0 for v in slab_checks.values()),
                            surface_pattern_geometry=dict(kind=graph['surface_pattern'],strut_diameter_mm=graph['surface_strut_diameter_mm'],
                                nominal_face_count=len(graph['surface_faces']),face_probes=surface_checks,
                                supported_faces_per_slab=surface_slabs,
                                every_slab_has_a_supported_surface_face=all(n>0 for n in surface_slabs.values())),
                            station_checks=station_checks,proximal_contact_checks=contacts,local_M3_hole_checks=holes,
                            distal_cap_normal_capture=capture,
                            active_core_guard_scope='R1.61 along retained core minus full routed R0.86 lumen, outside explicit terminal cylinders; nominal 0.75mm guard',
                            old_missing_core_wall_guard_mm3=float(old_missing.volume()),new_missing_core_wall_guard_mm3=float(new_missing.volume()),
                            newly_missing_core_wall_guard_mm3=float((new_missing-old_missing).volume()),
                            stabilized_same_0_75mm_guard=dict(outer_radius_mm=1.612,full_lumen_exclusion_radius_mm=.862,
                                radial_shift_mm=.002,missing_material_mm3=float(stable_missing.volume()),
                                interpretation='same 0.75mm ring shifted outward 0.002mm to avoid coincident float32 lumen faces; not print tolerance'),
                            parent_TPU_volume_mm3=float(old.volume()),new_TPU_volume_mm3=float(S[key].volume()),
                            surface_protective_pattern_verified=False,physical_shortening_mm=None,physical_force_N=None))
    changed=[key for key,r in generation['sources'].items() if r['sha256']!=parent['sources'][key]['sha256'] or r['translation_source_mm']!=parent['sources'][key]['translation_source_mm']]
    frames_ok=all(c['centre_error_mm']<1e-8 and c['orthogonality_error']<1e-10 and abs(c['determinant']-1.)<1e-10 and c['contraction_Y_dot_path_tangent']>1-1e-10
                  for m in muscles for c in m['station_checks'])
    report=dict(scope=__doc__,generation_source=source.relative_to(A.ROOT).as_posix(),generation_sha256=A.hashed(source),
                script_sha256=A.hashed(Path(__file__).resolve()),fixation_reference_sha256=A.hashed(fixation_file),
                parent_sha256=A.hashed(parent_file),n_instances=len(S),sources=generation['sources'],mesh_census=census,
                changed_keys=changed,pair_count_checked=len(S)*(len(S)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
                muscles=muscles,paper_Y_frame_checks_pass=frames_ok,
                all_axial_slabs_have_supported_cells=all(m['all_axial_slabs_have_supported_cells'] for m in muscles),
                no_rigid_overlap_over_threshold=not any(v['volume_mm3']>.001 and not v['compliant'] for v in overlaps),
                no_wire_overlap_over_threshold=not any(w['hits_over_threshold'] for w in wires),
                free_terminal_constraints_required=parent['free_terminal_constraints_required'],intended_body_owners=parent['intended_body_owners'],
                GH_nut_keepers=parent['GH_nut_keepers'],GH_assembly_entries=parent['GH_assembly_entries'],
                GH_retention_revision=parent['GH_retention_revision'],GH_lip_relief=parent['GH_lip_relief'],
                stage_pass=False,actual_drive_range_deg=None,
                not_verified=['nonlinear deformed geometry, channel stiffness/friction and free-terminal/tendon equilibrium',
                              'printed mesh/surface covering, process, real clamp/clip strength and measured muscle behaviour',
                              'all arm motion/K1-K7 and physical actuation'])
    (DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Changed',changed,'frames',frames_ok,'rigid clear',report['no_rigid_overlap_over_threshold'],'wire clear',report['no_wire_overlap_over_threshold'],flush=True)
    for m in muscles:print(m['module'],'contacts',[(x['part'],x['new_0_05mm_probe_overlap_mm3']) for x in m['proximal_contact_checks']],
                           'capture',m['distal_cap_normal_capture'],'wall new/old',m['new_missing_core_wall_guard_mm3'],m['old_missing_core_wall_guard_mm3'],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
