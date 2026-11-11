"""Assemble thoracic anatomy and the current shoulder/elbow bone references."""
import csv
import shutil

import numpy as np

from build_s1_teres_origin_mount import ROOT, mesh_geom, read, save, sha
from build_channel_first_inventory import provenance


def main():
    concept = ROOT/'02_人体骨骼晶格肌肉方案'
    folder = concept/'输出/channel_first'
    support_folder = concept/'输入/骨骼/支撑参考'
    support_file = support_folder/'sources.json'
    policy_file = ROOT/'design/CHANNEL_FIRST_SUPPORT_POLICY.json'
    policy = read(policy_file)
    assert policy['grounded_structure'] == 'thorax_on_tabletop_base'
    assembly_file = concept/'输出/integration_review/s1_elbow_shared_humerus/assembly_check.json'
    old = read(assembly_file)['sources']
    selected_groups = {'sternum', 'thoracic_spine', 'ribs', 'costal_cartilage'}
    models = list(read(support_file)['models'])
    database = ROOT.parent/'isa_BP3D_4.0_obj_99'
    table_file = database/'isa_element_parts.txt'
    additions = []
    for row in csv.DictReader(table_file.open(encoding='utf-8-sig'), delimiter='\t'):
        name = row['name'].lower()
        left_part = name.startswith('left ') and (name.endswith(' rib') or name.endswith(' costal cartilage'))
        right_lower = name in ('right tenth rib', 'right eleventh rib', 'right twelfth rib')
        if not (left_part or right_lower):
            continue
        part_id = row['element file id']
        source = database/f'{part_id}.obj'
        with source.open(encoding='utf-8') as stream:
            header = ''.join(next(stream) for _ in range(15))
        assert f'# File ID : {part_id}' in header
        assert f'# English name : {name}'.lower() in header.lower()
        local_folder = concept/'输入/骨骼/胸廓补充'
        local_folder.mkdir(exist_ok=True)
        destination = local_folder/source.name
        if destination.exists():
            assert sha(destination) == sha(source)
        else:
            shutil.copyfile(source, destination)
        item = dict(id=part_id, label=name, categories=['costal_cartilage' if 'cartilage' in name else 'ribs'],
                    sha256=sha(destination), local_file=destination.relative_to(ROOT).as_posix())
        models.append(item)
        additions.append(dict(**item, external_source=str(source)))
    save(concept/'输入/骨骼/胸廓补充/sources.json', dict(
        source_table=provenance(table_file) if table_file.is_relative_to(ROOT) else
        dict(file=str(table_file), sha256=sha(table_file)), models=additions))
    sources = {}
    regions = {}
    all_vertices = []
    for row in models:
        if not selected_groups.intersection(row['categories']):
            continue
        path = ROOT/row['local_file'] if 'local_file' in row else support_folder/f"{row['id']}.obj"
        assert sha(path) == row['sha256']
        vertices, faces = mesh_geom.load_obj(str(path))
        assert np.isfinite(vertices).all() and len(faces) > 0
        key = f"thorax_reference:{row['id']}"
        sources[key] = dict(file=path.relative_to(ROOT).as_posix(), sha256=sha(path),
                            translation_source_mm=[0., 0., 0.])
        regions[key] = dict(name=row['label'], categories=row['categories'],
                            bounds_source_mm=[*vertices.min(axis=0), *vertices.max(axis=0)],
                            vertex_count=len(vertices), triangle_count=len(faces))
        all_vertices.append(vertices)
    for key in ('bone:clavicle', 'bone:scapula', 'bone:humerus', 'bone:forearm'):
        assert sha(ROOT/old[key]['file']) == old[key]['sha256']
        sources[key] = old[key]
    vertices = np.vstack(all_vertices)
    bounds = [*vertices.min(axis=0), *vertices.max(axis=0)]
    save(folder/'thorax_arm_reference.json', dict(
        scope='Thorax anatomy and shoulder/elbow bone placement reference, without historical bench fixture.',
        sources=sources, n_instances=len(sources), support_inventory=provenance(support_file),
        arm_bone_inventory=provenance(assembly_file), support_policy=provenance(policy_file),
        coordinate_frame='Original BodyParts3D mm; no scaling, mirroring, or registration adjustment',
        thorax_regions=regions, thorax_bounds_source_mm=bounds,
        supplemental_thorax_inventory=provenance(concept/'输入/骨骼/胸廓补充/sources.json'),
        required_but_missing_in_input=[],
        structure_to_design=['thoracic spine-to-tabletop base', 'rib/cartilage/sternum fixed assembly',
                             'SC and scapular support', 'motor and conduit reaction mounts'],
        reference_only=True, assembled_mechanical_support_verified=False,
        collision_validated=False, drive_validated=False))
    print('thorax parts', len(regions), 'total scene', len(sources), 'thorax bounds', bounds)


if __name__ == '__main__':
    main()
