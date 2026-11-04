"""Export the corrected shoulder-elbow candidate for the static browser viewer."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
import mesh_geom
from audit_s1_extensor_guide_body import load_relocated_assembly
from build_s1_teres_origin_mount import A, load_component

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def read(path):return json.loads(path.read_text(encoding='utf8'))

def category(key):
    if key.startswith('bone:'):return 'bones'
    if key.startswith('tpu:'):return 'muscles'
    if key.startswith(('tendon:','wire:')) or key.endswith('/tendon_wire'):return 'wires'
    return 'hardware'

def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    assembly_file=base/'extensor_late_reaction_audit.json'
    bridge_file=base/'elbow_relocated_bridge_screen.json'
    audit=read(assembly_file)
    if sha(bridge_file)!=audit['bridge_sha256']:raise ValueError('Assembly report changed')
    joint=next(row['combined_key_and_seat'] for row in audit['planar_bearing_trials']
               if row.get('combined_key_and_seat', {}).get('sources'))
    pin=joint['retention_and_tool']['transverse_pin_trials']['selected']
    guide_file=base/'extensor_guide_body_audit.json'
    guide_audit=read(guide_file)
    channels=read(base/'elbow_relocated_channel_cores/channel_core_check.json')
    parts_by_key=load_relocated_assembly(
        base, read(bridge_file), channels, joint['sources']['tpu'],
        pin['axial_fit']['blind_thread_trial']['cap_source'])
    stop=next(row['source'] for row in guide_audit['detachable_guide']['assembly_readback']
              ['washer_tongue_trials'] if 'source' in row)
    for key, source in {
        'extensor:dual_key_guide':joint['sources']['guide'],
        'extensor:late_reaction_seat':joint['sources']['seat'],
        'extensor:retention_pin':pin['sources']['pin'],
        'elbow:FJ1477_proximal_nut_washer_right':stop,
    }.items():
        parts_by_key[key]=load_component(source)
    for channel in channels['channels']:
        label=channel['label']
        start,end=channels['choice'][f'{label}_endpoints_source_mm']
        middle=[guide_audit['guide_source_mm']] if label=='extensor' else channel['path_source_mm']
        parts_by_key[f'wire:elbow_{label}']=A.D.outer_envelope(np.array([start,*middle,end]),.3)
    source_parts=[]
    for key,solid in parts_by_key.items():
        vertices,faces=mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
        source_parts.append((key,vertices,faces))
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
    mesh_file=destination/'scene.bin';mesh_file.write_bytes(binary)
    manifest=dict(schema_version=1,candidate='S1 shoulder-elbow mechanical candidate',scope='Static neutral pose with nominal elbow cable paths. Motion and driven performance are outside this display.',
        assembly_source=assembly_file.relative_to(ROOT).as_posix(),assembly_sha256=sha(assembly_file),
        bridge_sha256=sha(bridge_file),guide_audit_sha256=sha(guide_file),
        exporter_sha256=sha(Path(__file__)),mesh_file='scene.bin',mesh_sha256=sha(mesh_file),mesh_bytes=len(binary),units='m',
        source_origin_mm=pivot.tolist(),bounds_m=[((lower-pivot)*.001).tolist(),((upper-pivot)*.001).tolist()],
        part_count=len(parts),triangle_count=sum(p['triangle_count'] for p in parts),
        decimated=False,maximum_coordinate_roundtrip_error_mm=maximum_error,
        category_counts={name:sum(p['category']==name for p in parts) for name in ['bones','muscles','wires','hardware']},parts=parts)
    (destination/'scene.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:manifest[k] for k in ['part_count','triangle_count','mesh_bytes','maximum_coordinate_roundtrip_error_mm','category_counts']},indent=2))

if __name__=='__main__':main()
