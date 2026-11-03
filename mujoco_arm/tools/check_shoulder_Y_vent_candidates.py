"""Independent readback and protected-stock checks of finite vent trials.

No selection or adoption is performed. Contact, material response, nozzle,
bridging, strength and muscle driving remain unverified.
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

DEST=A.DESIGN/'s1_y_muscle_pair/vent_experiment'

def run():
    source=DEST/'vent_trials.json';trial=json.loads(source.read_text(encoding='utf8'))
    assembly=json.loads((A.ROOT/trial['assembly_source']).read_text(encoding='utf8'))
    if trial['assembly_sha256']!=A.hashed(A.ROOT/trial['assembly_source']):raise ValueError('Stale assembly')
    for rel,key in [('build_shoulder_Y_vent_candidates.py','script_sha256'),('build_shoulder_Y_print_inspection.py','export_helper_sha256')]:
        if trial[key]!=A.hashed(A.ROOT/'mujoco_arm/tools'/rel):raise ValueError('Changed dependency')
    routefile=A.ROOT/trial['routes_source']
    if trial['routes_sha256']!=A.hashed(routefile):raise ValueError('Changed routes')
    routes=json.loads(routefile.read_text(encoding='utf8'))['routes']
    solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0]
        for k,r in assembly['sources'].items()}
    reference= A.DESIGN/'s1_GH_assembly_review/compact_recessed_keepers/assembly_check.json'
    parent=json.loads(reference.read_text(encoding='utf8'))
    void_roles=[];protected_by_module={}
    for module in ['infraspinatus','subscapularis']:
        key='tpu:'+module;full=solids[key];r=next(r for r in routes if r['module']==module)
        paths=[r['core_path_source_mm']]+[r['connectors'][s]['path_source_mm'] for s in ['origin','distal']]
        tubes={radius:A.I._manifold().Manifold.batch_boolean([A.D.outer_envelope(np.array(p),radius) for p in paths],A.I._manifold().OpType.Add) for radius in [1.8,3.3]}
        end_pieces=[]
        for port in r['ports'].values():
            c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source'])
            end_pieces.append(rod(c-n*8.,c+n*3.,12.))
        ends=A.I._manifold().Manifold.batch_boolean(end_pieces,A.I._manifold().OpType.Add)
        protected_by_module[module]=(ends+tubes[1.8])^full
        prior_row=parent['sources'][key]
        prior=A.load_checked(A.ROOT/prior_row['file'],prior_row['translation_source_mm'],prior_row['sha256'])[0]
        prior_voids=[positive_void(s) for s in prior.decompose() if s.volume() < -1e-6]
        prior_empty=A.I._manifold().Manifold.batch_boolean(prior_voids,A.I._manifold().OpType.Add)
        for shell in full.decompose():
            if shell.volume()>=-1e-6:continue
            void=positive_void(shell);volume=void.volume()
            void_roles.append(dict(module=module,volume_mm3=float(volume),bounds_source_mm=list(shell.bounding_box()),
                fraction_in_R1_8_channel_stock_region=A.I.intersection_volume(void,tubes[1.8])/volume,
                fraction_in_R3_3_retained_transition_region=A.I.intersection_volume(void,tubes[3.3])/volume,
                fraction_in_end_protection=A.I.intersection_volume(void,ends)/volume,
                overlap_with_prior_bone_baseline_void_fraction=A.I.intersection_volume(void,prior_empty)/volume,
                interpretation='geometric provenance only; not proof of functional intent or printability'))
    results=[]
    for row in trial['results']:
        key='tpu:'+row['module'];original=solids[key];mesh=row['mesh']
        actual,_,meta=A.load_checked(A.ROOT/mesh['file'],mesh['translation_source_mm'],mesh['sha256'])
        positive=[float(p.volume()) for p in actual.decompose() if p.volume()>1e-6]
        negative=[float(p.volume()) for p in actual.decompose() if p.volume() < -1e-6]
        hits=[]
        for other,s in solids.items():
            if other==key or is_compliant(other):continue
            value=A.I.intersection_volume(actual,s)
            if value>.001:hits.append(dict(part=other,volume_mm3=value))
        wire_hits=[]
        for r in routes:
            wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)
            value=A.I.intersection_volume(actual,wire)
            if value>.001:wire_hits.append(dict(head=r['head'],volume_mm3=value))
        protection=protected_by_module[row['module']]
        protected_loss=float((protection-actual).volume())
        aperture_checks=[]
        for step in row['steps']:
            if 'selected' not in step:continue
            chosen=step['selected'];p,q=np.array(chosen['path_source_mm'])
            # Finite cylinder replay in the unchanged source coordinate frame.
            aperture=rod(p,q,chosen['radius_mm'])
            aperture_checks.append(dict(radius_mm=chosen['radius_mm'],readback_material_remaining_in_cut_mm3=A.I.intersection_volume(aperture,actual)))
        geom_ok=(not hits and not wire_hits and meta['boundary_edges']==meta['nonmanifold_edges']==0 and len(positive)==1
                 and not negative and protected_loss<=.001 and all(c['readback_material_remaining_in_cut_mm3']<=.001 for c in aperture_checks))
        out=dict(module=row['module'],aperture_diameter_candidate_mm=row['aperture_diameter_candidate_mm'],
            positive_component_volumes_mm3=positive,enclosed_void_volumes_mm3=negative,
            protected_stock_loss_signed_mm3=protected_loss,report_volume_threshold_mm3=.001,
            readback_added_material_signed_mm3=float((actual-original).volume()),
            readback_removed_material_signed_mm3=float((original-actual).volume()),
            signed_volume_note='tiny negative Boolean volumes are numerical residuals, not physical negative material',
            rigid_hits_over_threshold=hits,wire_hits_over_threshold=wire_hits,aperture_replays=aperture_checks,
            rigid_parts_checked=sum(not is_compliant(k) for k in solids if k!=key),active_wire_paths_checked=len(routes),
            geometry_opening_checks_pass=geom_ok,
            rejected_reason=None if geom_ok else 'finite opening search failed or independent checks did not establish complete opening',
            adopted_as_current=False,manufacturing_or_drive_verified=False)
        results.append(out)
        print(row['module'],row['aperture_diameter_candidate_mm'],'geometry',geom_ok,'protected loss',protected_loss,flush=True)
    output=dict(scope=__doc__,trial_report_sha256=A.hashed(source),assembly_sha256=trial['assembly_sha256'],
        script_sha256=A.hashed(Path(__file__).resolve()),bone_reference_sha256=A.hashed(reference),
        source_void_role_probes=void_roles,results=results,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'independent_checks.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
