"""Bring original elbow interface parts into the new bones without moving muscles."""
from itertools import combinations
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'unified_SC_AC/assembly.json';parent=read(parent_file)
    original_file=ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus/assembly_check.json'
    original=read(original_file)
    baseline_file=original_file.parent/'anatomical_elbow_baseline.json';baseline=read(baseline_file)
    assert baseline['original_assembly_sha256']==provenance(original_file)['sha256']
    ownership_file=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/assembly.json'
    ownership={r['key']:r for r in read(ownership_file)['instances']}
    existing={**parent['sources'],**parent['joint_context_sources']}
    added={k:v for k,v in original['sources'].items() if k.startswith('elbow:') and k not in existing}
    solids={k:load_component(v) for k,v in {**existing,**added}.items()}
    hits=[];count=0
    for a,b in combinations(solids,2):
        if a not in added and b not in added:continue
        count+=1;volume=float(A.I.intersection_volume(solids[a],solids[b]))
        if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
    owners={k:dict(legacy_owner=ownership.get(k.removeprefix('elbow:'),{}).get('owner'),
                   source=provenance(ownership_file)) for k in added}
    display={f'tpu:elbow_{label}':row['display_source'] for label,row in baseline['muscles'].items()}
    muscle_records={label:dict(anatomical_id=row['anatomical_id'],
        pad_centres_source_mm=row['pad_centres_source_mm'],wire_path_source_mm=row['wire_path_source_mm'],
        source_definition=dict(file=row['definition'],sha256=row['definition_sha256']),
        status='original_geometry_preserved_no_moving_attachment_or_drive_approval')
        for label,row in baseline['muscles'].items()}
    dest=base/'elbow_interface_integration';dest.mkdir(exist_ok=True)
    save(dest/'assembly.json',dict(sources={**existing,**added},added_interface_sources=added,
        display_only_TPU_sources=display,legacy_owners=owners,muscles=muscle_records,
        parent=provenance(parent_file),original_rigid_assembly=provenance(original_file),
        original_muscle_baseline=provenance(baseline_file),ownership=provenance(ownership_file),
        generator=provenance(ROOT/'tools/restore_channel_first_elbow_interfaces.py'),
        changed_rigid_pairs_checked=count,neutral_interferences=hits,
        status='original_elbow_interfaces_restored_fit_pending' if hits else 'neutral_interface_geometry_candidate',
        scope='Original caps, bone frames, fasteners, keepers and swivel packs restored at exact original coordinates. Muscles and wire routes unchanged. Display TPU may not be a valid analysis solid. No cam, relocated muscle or prescribed swivel coupling adopted. Passive cap pivots and fixation need explicit modelling.'))
    save(dest/'view_scene.json',dict(sources={**existing,**added,**display},
        status='original_elbow_interface_context_not_motion_approved'))
    print('Restored',len(added),'rigid interfaces;',count,'changed pairs; hits',hits,flush=True)


if __name__=='__main__':main()
