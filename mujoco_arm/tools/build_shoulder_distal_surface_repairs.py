"""Small explicit TPU junction-stock candidates at measured PLC crossings.

Only add recorded local spheres to the complete muscle; do not move bones,
anchors or channels, delete hooks, fill arbitrary shells, or replace S1-Y2.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from build_shoulder_Y_print_inspection import write_readback
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_tpu_distal_surface_revision'

def run():
    selection_file=A.ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=json.loads(selection_file.read_text(encoding='utf8'))
    parent_file=A.ROOT/selection['assembly_report'];parent=json.loads(parent_file.read_text(encoding='utf8'))
    if A.hashed(parent_file)!=selection['assembly_sha256']:raise ValueError('Stale selected model')
    diagnostic_file=parent_file.parent/'volume_mesh/diagnostics.json';diagnostic=json.loads(diagnostic_file.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in parent['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=json.loads(route_file.read_text(encoding='utf8'))['routes']
    DEST.mkdir(parents=True,exist_ok=True);attempts=[]
    for failure in diagnostic['localized_failures']:
        module=failure['module'];key='tpu:'+module;source=parent['sources'][key];old=S[key];origin=np.array(source['translation_source_mm'])
        points=[]
        for crossing in failure['strict_nonadjacent_segment_facet_crossings']:
            p=np.array(crossing['intersection']['point_local_mm'])+origin
            if not any(np.linalg.norm(p-q)<.001 for q in points):points.append(p)
        own_route=next(r for r in routes if r['module']==module)
        lumen=A.D.outer_envelope(np.array(own_route['full_wire_path_source_mm']),.86)
        for radius in [.02,.05,.1,.2]:
            folder=DEST/module/('radius_'+f'{radius:g}');folder.mkdir(parents=True,exist_ok=True)
            patches=[A.I._manifold().Manifold.sphere(radius,24).translate(p) for p in points]
            stock=A.I._manifold().Manifold.batch_boolean(patches,A.I._manifold().OpType.Add)
            candidate=old+stock
            restored,mesh=write_readback(candidate.translate(-origin),folder/'TPU_distal_surface_candidate.stl')
            actual=restored.translate(origin);mesh['translation_source_mm']=origin.tolist()
            raw,faces=mesh_geom.read_stl(A.ROOT/mesh['file']);unique,inverse=np.unique(raw,axis=0,return_inverse=True)
            exact=mesh_geom.edge_census(inverse[faces]);V,F=mesh_geom.weld(raw,faces);imported=mesh_geom.edge_census(F)
            pieces=actual.decompose();positive=[float(p.volume()) for p in pieces if p.volume()>1e-6];negative=[float(p.volume()) for p in pieces if p.volume() < -1e-6]
            hits=[];soft=[]
            for other,body in S.items():
                if other==key:continue
                volume=A.I.intersection_volume(actual,body)
                if volume>.001:(soft if is_compliant(other) else hits).append(dict(part=other,volume_mm3=volume))
            wire_hits=[]
            for r in routes:
                wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3);volume=A.I.intersection_volume(actual,wire)
                if volume>.001:wire_hits.append(dict(head=r['head'],volume_mm3=volume))
            row=dict(module=module,patch_radius_mm=radius,patch_centres_source_mm=[p.tolist() for p in points],
                source=source,mesh=mesh,exact_boundary_edges=exact[1],exact_nonmanifold_edges=exact[2],
                imported_boundary_edges=imported[1],imported_nonmanifold_edges=imported[2],
                positive_component_volumes_mm3=positive,enclosed_void_volumes_mm3=negative,
                added_material_signed_mm3=float((actual-old).volume()),removed_original_material_signed_mm3=float((old-actual).volume()),
                full_routed_lumen_overlap_mm3=A.I.intersection_volume(actual,lumen),rigid_hits_over_threshold=hits,
                undeformed_TPU_hits_over_threshold=soft,active_wire_hits_over_threshold=wire_hits,
                static_preflight_pass=len(positive)==1 and not negative and not hits and not wire_hits and imported[1]==imported[2]==0,
                manufacturing_or_flexible_performance_verified=False,adopted_as_current=False,stage_pass=False)
            attempts.append(row)
            assembly=dict(scope='One complete-muscle local junction candidate; only its changed/context pairs checked, other 190 sources unchanged.',
                sources=dict(parent['sources']),parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
                source_change=row,stage_pass=False,actual_drive_range_deg=None)
            assembly['sources'][key]=mesh
            (folder/'candidate_assembly.json').write_text(json.dumps(assembly,indent=2)+'\n',encoding='utf8')
            print(module,radius,'added',row['added_material_signed_mm3'],'exact nm',exact[2],'static',row['static_preflight_pass'],flush=True)
    report=dict(scope=__doc__,parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
        diagnostic_source=diagnostic_file.relative_to(A.ROOT).as_posix(),diagnostic_sha256=A.hashed(diagnostic_file),
        script_sha256=A.hashed(Path(__file__).resolve()),export_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_Y_print_inspection.py'),
        attempts=attempts,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'junction_repair_trials.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
