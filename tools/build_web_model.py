"""Export selected S1 geometry for a static browser viewer, without decimation."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
import mesh_geom

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def read(path):return json.loads(path.read_text(encoding='utf8'))

def category(key):
    if key.startswith('bone:'):return 'bones'
    if key.startswith('tpu:'):return 'muscles'
    if key.startswith(('tendon:','wire:')) or key.endswith('/tendon_wire'):return 'wires'
    return 'hardware'

def main():
    pointer=ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=read(pointer)
    assembly_file=ROOT/selection['assembly_report'];assembly=read(assembly_file)
    if sha(assembly_file)!=selection['assembly_sha256']:raise ValueError('Selected report changed')
    source_parts=[]
    for key,row in assembly['sources'].items():
        path=ROOT/row['file']
        if sha(path)!=row['sha256']:raise ValueError('Source hash mismatch: '+key)
        vertices,faces=mesh_geom.load_obj(str(path)) if path.suffix.lower()=='.obj' else mesh_geom.read_stl(path)
        unique,inverse=np.unique(vertices,axis=0,return_inverse=True)
        faces=inverse[faces]
        source_parts.append((key,row,unique+np.array(row['translation_source_mm'] or [0,0,0]),faces))
    lower=np.min([v.min(axis=0) for _,_,v,_ in source_parts],axis=0)
    upper=np.max([v.max(axis=0) for _,_,v,_ in source_parts],axis=0)
    pivot=(lower+upper)/2
    destination=ROOT/'docs/viewer';destination.mkdir(parents=True,exist_ok=True)
    binary=bytearray();parts=[];maximum_error=0.
    for key,row,vertices,faces in source_parts:
        expected=(vertices-pivot)*.001;positions=np.asarray(expected,dtype='<f4')
        error=float(np.max(np.abs(positions.astype(float)-expected)))*1000;maximum_error=max(maximum_error,error)
        indices=np.asarray(faces,dtype='<u2' if len(vertices)<=65535 else '<u4')
        position_offset=len(binary);binary.extend(positions.tobytes())
        index_offset=len(binary);binary.extend(indices.tobytes())
        binary.extend(b'\0'*((-len(binary))%4))
        parts.append(dict(key=key,category=category(key),source_file=Path(row['file']).as_posix(),source_sha256=row['sha256'],
            vertex_count=len(vertices),triangle_count=len(faces),position_byte_offset=position_offset,
            index_byte_offset=index_offset,index_component_bytes=indices.dtype.itemsize,
            maximum_coordinate_roundtrip_error_mm=error))
    mesh_file=destination/'scene.bin';mesh_file.write_bytes(binary)
    manifest=dict(schema_version=1,candidate='S1-Y2',scope='Static display geometry, not a motion, contact, strength or TPU deformation simulation.',
        assembly_source=Path(selection['assembly_report']).as_posix(),assembly_sha256=sha(assembly_file),selection_sha256=sha(pointer),
        exporter_sha256=sha(Path(__file__)),mesh_file='scene.bin',mesh_sha256=sha(mesh_file),mesh_bytes=len(binary),units='m',
        source_origin_mm=pivot.tolist(),bounds_m=[((lower-pivot)*.001).tolist(),((upper-pivot)*.001).tolist()],
        part_count=len(parts),triangle_count=sum(p['triangle_count'] for p in parts),
        source_triangles_preserved=True,decimated=False,maximum_coordinate_roundtrip_error_mm=maximum_error,
        category_counts={name:sum(p['category']==name for p in parts) for name in ['bones','muscles','wires','hardware']},parts=parts)
    (destination/'scene.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:manifest[k] for k in ['part_count','triangle_count','mesh_bytes','maximum_coordinate_roundtrip_error_mm','category_counts']},indent=2))

if __name__=='__main__':main()
