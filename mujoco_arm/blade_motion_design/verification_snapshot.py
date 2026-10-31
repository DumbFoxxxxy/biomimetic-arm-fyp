"""Bind audit results to the exact XML, asset bytes and target trajectories."""
import hashlib,xml.etree.ElementTree as ET
from pathlib import Path
from arm_targets import target_actions

def fingerprint(root):
    xml_hash=hashlib.sha256(ET.tostring(root,encoding='utf-8')).hexdigest();assets={}
    for mesh in root.findall('./asset/mesh'):
        path=Path(mesh.get('file'))
        if not path.is_absolute():path=Path(root.find('compiler').get('meshdir','.'))/path
        assets[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    trajectory=hashlib.sha256()
    for name,points in target_actions().items():trajectory.update(name.encode());trajectory.update(points.tobytes())
    return {'xml_sha256':xml_hash,'mesh_sha256':assets,'trajectories_sha256':trajectory.hexdigest()}
