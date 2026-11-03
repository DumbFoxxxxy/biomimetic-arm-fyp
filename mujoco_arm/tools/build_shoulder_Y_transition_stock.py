"""Explicit design alternatives for classified closed transition pockets.

Infill is a recorded material change, not export cleanup. Only closed voids
inside the R3.3 retained transition region and outside routed lumens/end
interfaces qualify; collisions and readback must be checked independently.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_Y_vent_candidates import positive_void,negative_volume
from build_shoulder_Y_print_inspection import write_readback
from build_shoulder_shaft_collar import rod

DEST=A.DESIGN/'s1_y_muscle_pair/transition_stock_comparison'

def run():
    assembly_file=A.DESIGN/'s1_y_muscle_pair/assembly_check.json'
    assembly=json.loads(assembly_file.read_text(encoding='utf8'))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(routefile.read_text(encoding='utf8'))['routes']
    outward_file=A.DESIGN/'s1_y_muscle_pair/vent_experiment/outward_start_revision/outward_trials.json'
    outward=json.loads(outward_file.read_text(encoding='utf8'))
    if outward['assembly_sha256']!=A.hashed(assembly_file):raise ValueError('Stale outward source')
    DEST.mkdir(parents=True,exist_ok=True);results=[]
    for module in ['infraspinatus','subscapularis']:
        row=assembly['sources']['tpu:'+module];original,_,_=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
        route=next(r for r in routes if r['module']==module)
        paths=[route['core_path_source_mm']]+[route['connectors'][s]['path_source_mm'] for s in ['origin','distal']]
        transition=A.I._manifold().Manifold.batch_boolean([A.D.outer_envelope(np.array(p),3.3) for p in paths],A.I._manifold().OpType.Add)
        lumen=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.86)
        ends=[]
        for port in route['ports'].values():
            c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source']);ends.append(rod(c-n*8.,c+n*3.,12.))
        end=A.I._manifold().Manifold.batch_boolean(ends,A.I._manifold().OpType.Add)
        options=[('closed_pockets_infilled',original,None)]
        if module=='subscapularis':
            vent=next(r for r in outward['results'] if r['aperture_diameter_candidate_mm']==.8)
            meta=vent['mesh'];body,_,_=A.load_checked(A.ROOT/meta['file'],meta['translation_source_mm'],meta['sha256'])
            options.append(('bounded_0_8mm_vents_plus_residual_infill',body,meta))
        for name,body,precursor in options:
            changes=[];candidate=body
            for shell in body.decompose():
                if shell.volume()>=-1e-6:continue
                pocket=positive_void(shell);volume=float(pocket.volume())
                outside=float((pocket-transition).volume());lumen_hit=A.I.intersection_volume(pocket,lumen);end_hit=A.I.intersection_volume(pocket,end)
                if max(outside,lumen_hit,end_hit)>.000001:
                    raise ValueError('Pocket classification does not permit infill')
                # This Boolean is the deliberate design change. No other
                # shell is filled and no export repair runs afterwards.
                candidate=candidate+pocket
                changes.append(dict(empty_volume_infilled_mm3=volume,bounds_source_mm=list(shell.bounding_box()),
                    outside_R3_3_transition_region_mm3=outside,nominal_lumen_intersection_mm3=lumen_hit,
                    terminal_interface_intersection_mm3=end_hit,
                    rationale='closed non-route pocket in retained transition stock; intentionally replace void with TPU',
                    stiffness_and_deformation_effect_unverified=True))
            folder=DEST/module/name;folder.mkdir(parents=True,exist_ok=True)
            origin=np.array(row['translation_source_mm']);local,mesh=write_readback(candidate.translate(-origin),folder/'TPU_transition_trial.stl')
            actual=local.translate(origin);mesh['translation_source_mm']=origin.tolist()
            result=dict(module=module,variant=name,parent_TPU=row,precursor_mesh=precursor,explicit_infill_changes=changes,mesh=mesh,
                source_closed_void_volume_mm3=negative_volume(original),remaining_closed_void_volume_mm3=negative_volume(actual),
                added_material_signed_mm3=float((actual-original).volume()),removed_material_signed_mm3=float((original-actual).volume()),
                volume_change_percent=100.*float(actual.volume()-original.volume())/original.volume(),
                no_outer_material_removal_by_design=precursor is None,
                no_export_mesh_simplification=True,adopted_as_current=False,stage_pass=False)
            results.append(result)
            print(module,name,'void',result['remaining_closed_void_volume_mm3'],'added',result['added_material_signed_mm3'],
                  'removed',result['removed_material_signed_mm3'],flush=True)
    report=dict(scope=__doc__,assembly_source=assembly_file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(assembly_file),
        routes_source=routefile.relative_to(A.ROOT).as_posix(),routes_sha256=A.hashed(routefile),
        outward_source=outward_file.relative_to(A.ROOT).as_posix(),outward_sha256=A.hashed(outward_file),
        script_sha256=A.hashed(Path(__file__).resolve()),dependency_sha256={name:A.hashed(A.ROOT/'mujoco_arm/tools'/name)
            for name in ['build_shoulder_Y_vent_candidates.py','build_shoulder_Y_print_inspection.py']},
        results=results,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'transition_trials.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
