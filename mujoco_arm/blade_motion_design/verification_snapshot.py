"""Bind audit results to the exact XML, asset bytes and target trajectories."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import hashlib,xml.etree.ElementTree as ET
from project_runtime import WorkspacePath as Path
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
