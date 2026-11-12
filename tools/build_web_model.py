"""Export the shoulder, original elbow and integrated shoulder–elbow scenes."""
import argparse
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
import mesh_geom
from build_s1_teres_origin_mount import A

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def read(path):return json.loads(path.read_text(encoding='utf8'))

def category(key):
    if key.startswith(('bone:', 'thorax_reference:')):return 'bones'
    if key.startswith('tpu:'):return 'muscles'
    if key.startswith(('tendon:','wire:','wire_active:')) or key.endswith('/tendon_wire'):return 'wires'
    return 'hardware'

def load_parts(sources):
    parts = []
    for key, row in sources.items():
        path = ROOT/row['file']
        if sha(path) != row['sha256']:
            raise ValueError(f'Mesh source changed: {key}')
        vertices, faces = (mesh_geom.load_obj(str(path)) if path.suffix.lower() == '.obj'
                           else mesh_geom.read_stl(path))
        unique, inverse = np.unique(vertices, axis=0, return_inverse=True)
        parts.append((key, unique+np.array(row['translation_source_mm'] or [0, 0, 0]), inverse[faces]))
    return parts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--combined-only', action='store_true')
    args = parser.parse_args()
    integration = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review'
    if args.combined_only:
        write_combined(integration)
        return
    base = integration/'s1_elbow_shared_humerus'
    assembly = base/'assembly_check.json'
    baseline = read(base/'anatomical_elbow_baseline.json')
    if sha(assembly) != baseline['original_assembly_sha256']:
        raise ValueError('Original elbow source changed')
    elbow_sources = {key: source for key, source in read(assembly)['sources'].items()
                     if key.startswith(('elbow:', 'tpu:elbow_')) or key == 'bone:forearm'}
    preview = read(integration/'s1_elbow_preview/integration_audit.json')
    elbow_sources['bone:humerus'] = preview['original_elbow_humerus']
    parts = load_parts(elbow_sources)
    for label, row in baseline['muscles'].items():
        if elbow_sources[f'tpu:elbow_{label}'] != row['display_source']:
            raise ValueError(f'Elbow muscle source changed: {label}')
        cable = A.D.outer_envelope(np.array(row['wire_path_source_mm']), .3)
        vertices, faces = mesh_geom.manifold_mesh_to_arrays(cable.to_mesh64())
        parts.append((f'wire:elbow_{label}', vertices, faces))
    write_scene(parts, 'scene', 'Original anatomical elbow', assembly)

    selection = read(ROOT/'design/S1_CURRENT_CANDIDATE.json')
    shoulder = ROOT/selection['assembly_report']
    if sha(shoulder) != selection['assembly_sha256']:
        raise ValueError('Shoulder source changed')
    write_scene(load_parts(read(shoulder)['sources']), 'shoulder', 'S1-Y2 shoulder', shoulder)
    write_combined(integration)


def write_combined(integration):
    folder = integration / 's1_series_elbow'
    assembly_file = folder / 'assembly_FJ1504_local_interface_bidirectional.json'
    assembly = read(assembly_file)
    relief_file = folder / 'FJ1486_return_relief.json'
    relief = read(relief_file)
    muscles = dict(assembly['flexible_geometry_sources'])
    if muscles['tpu:elbow_flexor'] != relief['original']:
        raise ValueError('The flexor relief is based on a different muscle mesh')
    muscles['tpu:elbow_flexor'] = relief['selected']['mesh']
    if assembly['sources'].keys() & muscles.keys():
        raise ValueError('An integrated part is listed as both rigid and flexible')
    sources = {**assembly['sources'], **muscles}
    parts = load_parts(sources)
    audit_file = folder / 'timed_elbow_muscle_geometry.json'
    audit = read(audit_file)
    passive_file = ROOT / audit['passive']['file']
    if sha(passive_file) != audit['passive']['sha256']:
        raise ValueError('The elbow path reference changed after verification')
    passive = read(passive_file)
    elbow_file = ROOT / passive['assembly']['file']
    if sha(elbow_file) != passive['assembly']['sha256']:
        raise ValueError('The original elbow assembly changed')
    for label, row in read(elbow_file)['muscles'].items():
        cable = A.D.outer_envelope(np.array(row['wire_path_source_mm'])[1:-1], .3)
        vertices, faces = mesh_geom.manifold_mesh_to_arrays(cable.to_mesh64())
        parts.append((f'wire:elbow_{label}', vertices, faces))
    write_scene(parts, 'combined', 'Integrated S1 shoulder and series-slider elbow', assembly_file,
                scope='Static integrated neutral candidate with anatomical thorax reference and revised elbow flexor.'
                      ' Only local elbow wires are shown. Some shoulder terminals are pending integration.',
                evidence={'elbow_flexor_relief': relief_file, 'elbow_wire_paths': elbow_file})


def write_scene(source_parts, stem, candidate, assembly_file, *,
                scope='Separate static neutral assembly. Motion and driven performance are outside this display.',
                evidence=None):
    lower=np.min([v.min(axis=0) for _,v,_ in source_parts],axis=0)
    upper=np.max([v.max(axis=0) for _,v,_ in source_parts],axis=0)
    pivot=(lower+upper)/2
    destination=ROOT/'docs/viewer';destination.mkdir(parents=True,exist_ok=True)
    binary=bytearray();parts=[];maximum_error=0.
    for key,vertices,faces in source_parts:
        expected=(vertices-pivot)*.001;positions=np.asarray(expected,dtype='<f4')
        error=float(np.max(np.abs(positions.astype(float)-expected)))*1000;maximum_error=max(maximum_error,error)
        indices=np.asarray(faces,dtype='<u2' if len(vertices)<=65535 else '<u4')
        position_offset=len(binary);binary.extend(positions.tobytes())
        index_offset=len(binary);binary.extend(indices.tobytes())
        binary.extend(b'\0'*((-len(binary))%4))
        parts.append(dict(key=key,category=category(key),
            vertex_count=len(vertices),triangle_count=len(faces),position_byte_offset=position_offset,
            index_byte_offset=index_offset,index_component_bytes=indices.dtype.itemsize,
            maximum_coordinate_roundtrip_error_mm=error))
    mesh_file=destination/f'{stem}.bin';mesh_file.write_bytes(binary)
    manifest=dict(schema_version=1,candidate=candidate,scope=scope,
        assembly_source=assembly_file.relative_to(ROOT).as_posix(),assembly_sha256=sha(assembly_file),
        exporter_sha256=sha(Path(__file__)),mesh_file=mesh_file.name,mesh_sha256=sha(mesh_file),mesh_bytes=len(binary),units='m',
        source_origin_mm=pivot.tolist(),bounds_m=[((lower-pivot)*.001).tolist(),((upper-pivot)*.001).tolist()],
        part_count=len(parts),triangle_count=sum(p['triangle_count'] for p in parts),
        decimated=False,maximum_coordinate_roundtrip_error_mm=maximum_error,
        category_counts={name:sum(p['category']==name for p in parts) for name in ['bones','muscles','wires','hardware']},parts=parts)
    if evidence:
        manifest['evidence'] = {key: {'file': path.relative_to(ROOT).as_posix(), 'sha256': sha(path)}
                                for key, path in evidence.items()}
    (destination/f'{stem}.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:manifest[k] for k in ['part_count','triangle_count','mesh_bytes','maximum_coordinate_roundtrip_error_mm','category_counts']},indent=2))

if __name__=='__main__':main()
