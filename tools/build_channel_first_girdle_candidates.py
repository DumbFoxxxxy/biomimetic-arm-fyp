"""Extract shoulder-girdle endpoint candidates from existing anatomical patches."""
import argparse
import numpy as np

from build_s1_teres_origin_mount import ROOT, read, save
from build_channel_first_inventory import provenance


def representative(patch):
    points = patch['bone_points_source_mm']
    weights = patch['muscle_face_areas_mm2']
    assert len(points) > 0 and weights.sum() > 0
    mean = np.average(points, axis=0, weights=weights)
    return points[np.argmin(np.linalg.norm(points-mean, axis=1))].tolist()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upper-girdle', action='store_true')
    parser.add_argument('--arm-antagonists', action='store_true')
    args = parser.parse_args()
    if args.arm_antagonists and not args.upper_girdle:
        parser.error('--arm-antagonists requires --upper-girdle')
    concept = ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept'
    region_folder = concept/'attachment_regions'
    atlas_file = region_folder/'atlas_and_regions.json'
    atlas = {row['id']: row for row in read(atlas_file)['muscles']}
    choices = [('FJ1456', 'ribs', 1), ('FJ1459', 'ribs', 3),
               ('FJ1520', 'thoracic_spine', 1), ('FJ1536', 'thoracic_spine', 2)]
    if args.upper_girdle:
        choices.extend([('FJ1460', 'first_rib', 1), ('FJ1521', 'cervical_spine', 1),
                        ('FJ1532', 'cervical_spine', 1), ('FJ1537', 'cervical_spine', 1)])
    choices = [(head, origin, count, 'clavicle' if head in ('FJ1460', 'FJ1521') else 'scapula')
               for head, origin, count in choices]
    if args.arm_antagonists:
        choices.extend([('FJ1506', 'scapula', 1, 'humerus'), ('FJ1508', 'scapula', 1, 'humerus'),
                        ('FJ1447', 'clavicle', 1, 'humerus'), ('FJ1464', 'sternum', 2, 'humerus'),
                        ('FJ1488', 'scapula', 1, 'humerus'), ('FJ1478', 'scapula', 1, 'radius'),
                        ('FJ1479', 'scapula', 1, 'ulna')])
    routes = []
    inputs = []
    for head, origin_structure, divisions, insertion_structure in choices:
        row = atlas[head]
        patches = {}
        origin_body = origin_structure if origin_structure in ('clavicle', 'scapula', 'humerus') else 'thorax'
        insertion_body = 'forearm' if insertion_structure in ('radius', 'ulna') else insertion_structure
        for role, structure in (('origin', origin_structure), ('insertion', insertion_structure)):
            record = next(item for item in row['regions'] if item['role'] == role and item['structure'] == structure)
            path = region_folder/record['data']
            assert record['selected_face_count_2mm'] > 0
            with np.load(path) as data:
                patch = {key: data[key] for key in ('bone_points_source_mm', 'muscle_face_areas_mm2')}
            points = patch['bone_points_source_mm']
            assert len(points) == record['selected_face_count_2mm']
            # Subdivide a broad attachment patch by height. This is an
            # engineering sampling rule with no claim of muscle fiber mapping.
            order = np.argsort(points[:, 2])
            patches[role] = [
                {key: values[indices] for key, values in patch.items()}
                for indices in np.array_split(order, divisions)]
            inputs.append(provenance(path))
        for index, (proximal, distal) in enumerate(zip(patches['origin'], patches['insertion'])):
            routes.append(dict(
                id=f'{head}_g{index+1}', anatomical_id=head, region='girdle_candidate',
                module=row['name'], anatomical_origin_body=origin_body, distal_body=insertion_body,
                candidate_endpoints_source_mm=[representative(proximal), representative(distal)],
                candidate_origin_structure=origin_structure,
                status='unreviewed proximity-derived endpoint candidate',
                origin_patch_samples=len(proximal['bone_points_source_mm']),
                insertion_patch_samples=len(distal['bone_points_source_mm']),
                support='thorax/scapula attachment hardware and routed channel need design'))
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    filename = 'upper_girdle_endpoint_candidates.json' if args.upper_girdle else 'girdle_endpoint_candidates.json'
    if args.arm_antagonists:
        filename = 'arm_girdle_endpoint_candidates.json'
    save(folder/filename, dict(
        anatomical_patch_report=provenance(atlas_file), patch_inputs=inputs, routes=routes,
        included_muscles=[item[0] for item in choices],
        pending_upper_girdle=['mechanical cervical support tied to the grounded thorax',
                             'attachment hardware', 'routed channel definition'],
        scope='Candidate points sampled on previously projected bone surfaces. No clinical landmark, '
              'fiber direction, mechanical fixation, guide path or collision approval.'))
    print('girdle endpoint candidates', len(routes))


if __name__ == '__main__':
    main()
