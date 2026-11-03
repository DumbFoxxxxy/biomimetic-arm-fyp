"""Independent audit of bounded-start vents and explicit pocket infill.

Checks cover actual meshes, zero-pose rigid/wire fits, protected stock,
finite interface contacts, emitted lattice coverage and design replay.
They do not cover mechanical response, slicing, strength or driven motion.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_shaft_collar import rod
from build_shoulder_Y_vent_candidates import positive_void
from check_shoulder_front_liner_keepers import is_compliant

BASE=A.DESIGN/'s1_y_muscle_pair'
DEST=BASE/'transition_stock_comparison'

def load(row):return A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])[0]

def run():
    assembly_file=BASE/'assembly_check.json';assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    generation=json.loads((BASE/'candidate_generation.json').read_text(encoding='utf8'))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=json.loads(routefile.read_text(encoding='utf8'))['routes']
    outward_file=BASE/'vent_experiment/outward_start_revision/outward_trials.json'
    transition_file=DEST/'transition_trials.json'
    documents=[(outward_file,'build_shoulder_Y_outward_vents.py'),(transition_file,'build_shoulder_Y_transition_stock.py')]
    trials=[]
    for file,script in documents:
        d=json.loads(file.read_text(encoding='utf8'))
        if d['assembly_sha256']!=A.hashed(assembly_file) or d['routes_sha256']!=A.hashed(routefile):raise ValueError('Stale sources')
        if d['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools'/script):raise ValueError('Changed generator')
        if any(A.hashed(A.ROOT/'mujoco_arm/tools'/k)!=v for k,v in d['dependency_sha256'].items()):raise ValueError('Changed dependency')
        if 'outward_sha256' in d and d['outward_sha256']!=A.hashed(outward_file):raise ValueError('Stale precursor report')
        for r in d['results']:trials.append((file,r))
    solids={k:load(v) for k,v in assembly['sources'].items()};results=[]
    wire_solids=[(r['head'],A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)) for r in routes]
    for source,row in trials:
        module=row['module'];key='tpu:'+module;original=solids[key];actual=load(row['mesh'])
        route=next(r for r in routes if r['module']==module);parts=[]
        for p in route['ports'].values():
            c=np.array(p['inner_face_source_mm']);n=np.array(p['normal_source']);parts.append(rod(c-n*8.,c+n*3.,12.))
        for path in [route['core_path_source_mm']]+[route['connectors'][side]['path_source_mm'] for side in ['origin','distal']]:
            parts.append(A.D.outer_envelope(np.array(path),1.8))
        protected=original^A.I._manifold().Manifold.batch_boolean(parts,A.I._manifold().OpType.Add)
        if 'variant' in row:
            expected=original if row['precursor_mesh'] is None else load(row['precursor_mesh'])
            for shell in list(expected.decompose()):
                if shell.volume() < -1e-6:expected=expected+positive_void(shell)
            label=row['variant']
        else:
            expected=original
            for step in row['steps']:
                if 'selected' in step:
                    chosen=step['selected'];p,q=np.array(chosen['path_source_mm'])
                    expected=expected-rod(p,q,chosen['radius_mm'])
            label='bounded_start_vent_'+f'{row["aperture_diameter_candidate_mm"]:g}'
        replay_extra=float((actual-expected).volume());replay_lost=float((expected-actual).volume())
        protected_loss=float((protected-actual).volume());rigid=[];wires=[]
        for other,body in solids.items():
            if other==key or is_compliant(other):continue
            v=A.I.intersection_volume(actual,body)
            if v>.001:rigid.append(dict(part=other,volume_mm3=v))
        for head,wire in wire_solids:
            v=A.I.intersection_volume(actual,wire)
            if v>.001:wires.append(dict(head=head,volume_mm3=v))
        selected=next(r for r in generation['attempts'] if r['module']==module and r['selected'])
        graph=selected['transported_frames'];nodes={tuple(n['index']):np.array(n['source_mm']) for n in graph['nodes']}
        core_slabs={i:0 for i in range(len(graph['stations'])-1)};surface_slabs=dict(core_slabs)
        for cell in graph['cells']:
            c=np.array(cell['centre_source_mm']);count=0
            for index in cell['corner_indices']:
                q=nodes[tuple(index)];supported=False
                for f in [.35,.65]:
                    probe=A.I._manifold().Manifold.sphere(.15,16).translate(c+(q-c)*f)
                    if A.I.intersection_volume(probe,actual)>.8*probe.volume():supported=True;break
                count+=supported
            if count>=3:core_slabs[cell['index'][2]]+=1
        for face in graph['surface_faces']:
            c=np.array(face['centre_source_mm']);count=0
            for index in face['corner_indices']:
                q=nodes[tuple(index)];supported=False
                for f in [.4,.7]:
                    probe=A.I._manifold().Manifold.sphere(.1,16).translate(c+(q-c)*f)
                    if A.I.intersection_volume(probe,actual)>.8*probe.volume():supported=True;break
                count+=supported
            if count>=2:surface_slabs[face['axial_slab']]+=1
        head=route['head'];origin_n=np.array(route['ports']['origin']['normal_source']);contact=[]
        for part,n in [('cap:'+head+'/origin',-origin_n),('origin:'+head+'/front_support',origin_n)]:
            contact.append(dict(part=part,normal_0_05mm_overlap_mm3=A.I.intersection_volume(solids[part].translate(n*.05),actual)))
        n=np.array(route['ports']['distal']['normal_source']);cap=solids['free_terminal:'+head]
        for distance in [-.2,.25]:
            contact.append(dict(part='free_terminal:'+head,normal_distance_mm=distance,
                overlap_mm3=A.I.intersection_volume(cap.translate(n*distance),actual)))
        pieces=actual.decompose();positive=[float(p.volume()) for p in pieces if p.volume()>1e-6];negative=[float(p.volume()) for p in pieces if p.volume() < -1e-6]
        keep_coverage=all(core_slabs.values()) and all(surface_slabs.values())
        positive_contacts=all(c.get('normal_0_05mm_overlap_mm3',c.get('overlap_mm3',0))>.001 for c in contact)
        replay_ok=max(replay_extra,replay_lost)<=.001
        static_ok=(len(positive)==1 and not rigid and not wires and protected_loss<=.001 and replay_ok and keep_coverage and positive_contacts)
        result=dict(module=module,variant=label,source_report=source.relative_to(A.ROOT).as_posix(),mesh_sha256=row['mesh']['sha256'],
            positive_components=len(positive),closed_void_volumes_mm3=negative,
            expected_design_replay_extra_signed_mm3=replay_extra,expected_design_replay_loss_signed_mm3=replay_lost,
            protected_original_stock_loss_signed_mm3=protected_loss,
            rigid_parts_checked=177,rigid_hits_over_threshold=rigid,wire_hits_over_threshold=wires,
            supported_internal_cells_per_slab=core_slabs,supported_surface_faces_per_slab=surface_slabs,
            all_axial_slabs_supported=keep_coverage,finite_contact_probes=contact,finite_contacts_encounter_material=positive_contacts,
            static_geometry_checks_pass=static_ok,all_pockets_resolved=not negative,
            complete_pocket_resolution_geometry_pass=static_ok and not negative,
            signed_volume_note='small signed numerical Boolean residuals remain reported; threshold is 0.001 mm3, not a manufacturing tolerance',
            stage_pass=False,adopted_as_current=False,actual_drive_range_deg=None)
        results.append(result)
        print(module,label,'static',static_ok,'closed voids',negative,flush=True)
    output=dict(scope=__doc__,assembly_sha256=A.hashed(assembly_file),routes_sha256=A.hashed(routefile),
        trial_report_sha256={f.relative_to(A.ROOT).as_posix():A.hashed(f) for f,_ in documents},
        script_sha256=A.hashed(Path(__file__).resolve()),generation_sha256=A.hashed(BASE/'candidate_generation.json'),
        all_191_source_hashes_match=True,results=results,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'independent_checks.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
