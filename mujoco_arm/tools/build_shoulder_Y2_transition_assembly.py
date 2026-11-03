"""Assemble a research candidate with explicitly infilled transition pockets.

Only two TPU sources change. Classified pockets are outside all routed
lumens and terminal interface regions. No printing or actuation pass.
"""
import itertools
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import audit_shoulder_retained_joint_fit as A
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_y2_transition_stock'

def run():
    parent_file=A.DESIGN/'s1_y_muscle_pair/assembly_check.json';parent=json.loads(parent_file.read_text(encoding='utf8'))
    generation_file=A.DESIGN/'s1_y_muscle_pair/transition_stock_comparison/transition_trials.json'
    checked_file=generation_file.parent/'independent_checks.json'
    generation=json.loads(generation_file.read_text(encoding='utf8'));checked=json.loads(checked_file.read_text(encoding='utf8'))
    if generation['assembly_sha256']!=A.hashed(parent_file) or checked['assembly_sha256']!=A.hashed(parent_file):raise ValueError('Changed parent')
    if checked['trial_report_sha256'][generation_file.relative_to(A.ROOT).as_posix()]!=A.hashed(generation_file):raise ValueError('Stale trial audit')
    if checked['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/check_shoulder_Y_transition_stock.py'):raise ValueError('Changed audit')
    if generation['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_Y_transition_stock.py'):raise ValueError('Changed generator')
    sources=dict(parent['sources']);muscles=[]
    for row in generation['results']:
        if row['variant']!='closed_pockets_infilled':continue
        verified=next(r for r in checked['results'] if r['module']==row['module'] and r['variant']==row['variant'])
        if not verified['complete_pocket_resolution_geometry_pass'] or verified['mesh_sha256']!=row['mesh']['sha256']:raise ValueError('Candidate not verified')
        sources['tpu:'+row['module']]=row['mesh']
        muscles.append(dict(module=row['module'],variant=row['variant'],verification=verified,
            explicit_material_changes=row['explicit_infill_changes'],added_material_signed_mm3=row['added_material_signed_mm3'],
            removed_material_signed_mm3=row['removed_material_signed_mm3'],volume_change_percent=row['volume_change_percent'],
            physical_shortening_mm=None,physical_force_N=None))
    if len(muscles)!=2:raise ValueError('Need both complete muscle candidates')
    solids={};census={}
    for key,row in sources.items():
        solid,_,meta=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256']);solids[key]=solid
        pieces=solid.decompose();positive=[p for p in pieces if p.volume()>1e-6]
        if len(positive)!=1:raise ValueError('Disconnected assembly object '+key)
        census[key]=dict(boundary=meta['boundary_edges'],nonmanifold=meta['nonmanifold_edges'],positive_components=len(positive),
            enclosed_void_shell_volumes_mm3=[float(p.volume()) for p in pieces if p.volume() < -1e-6])
    overlaps=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=volume,compliant=is_compliant(a) or is_compliant(b)))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=json.loads(routefile.read_text(encoding='utf8'))['routes'];wires=[]
    for r in routes:
        wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)
        hits={k:A.I.intersection_volume(wire,s) for k,s in solids.items()}
        wires.append(dict(head=r['head'],volumes_mm3=hits,hits_over_threshold=[dict(part=k,volume_mm3=v) for k,v in hits.items() if v>.001]))
    changed=[k for k in sources if sources[k]['sha256']!=parent['sources'][k]['sha256'] or sources[k]['translation_source_mm']!=parent['sources'][k]['translation_source_mm']]
    if set(changed)!={'tpu:infraspinatus','tpu:subscapularis'}:raise ValueError('Unexpected source changes')
    rigid_ok=not any(r['volume_mm3']>.001 and not r['compliant'] for r in overlaps);wire_ok=not any(r['hits_over_threshold'] for r in wires)
    if not rigid_ok or not wire_ok:raise ValueError('Integrated static fit fails')
    report=dict(scope=__doc__,parent_report=parent_file.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(parent_file),
        generation_source=generation_file.relative_to(A.ROOT).as_posix(),generation_sha256=A.hashed(generation_file),
        independent_muscle_checks=checked_file.relative_to(A.ROOT).as_posix(),independent_muscle_checks_sha256=A.hashed(checked_file),
        script_sha256=A.hashed(Path(__file__).resolve()),n_instances=len(sources),sources=sources,mesh_census=census,
        changed_keys=changed,pair_count_checked=len(sources)*(len(sources)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
        muscles=muscles,no_rigid_overlap_over_threshold=rigid_ok,no_wire_overlap_over_threshold=wire_ok,
        all_axial_slabs_have_supported_cells=all(m['verification']['all_axial_slabs_supported'] for m in muscles),
        unchanged_bone_and_anchor_sources=True,frame_and_lattice_graph_reference=parent['generation_source'],
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],intended_body_owners=parent['intended_body_owners'],
        stage_pass=False,actual_drive_range_deg=None,
        not_verified=['TPU stiffness, nonlinear deformation, friction and free-terminal/tendon equilibrium',
            'slicer/device/material process, clamp strength and measured muscle behaviour','all arm K1-K7 or physical actuation'])
    for key in ['GH_nut_keepers','GH_assembly_entries','GH_retention_revision','GH_lip_relief']:report[key]=parent[key]
    DEST.mkdir(parents=True,exist_ok=True);(DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('S1-Y2',len(sources),'sources;',report['pair_count_checked'],'pairs; rigid/wire',rigid_ok,wire_ok,flush=True)

if __name__=='__main__':run()
