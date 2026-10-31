"""Use full existing all-pair motion audit with isolated M7 sources."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import sys,json
sys.argv.append('--muscle-mounts')
import audit_joint_model as audit
from project_runtime import WorkspacePath as Path
DEST=Path(__file__).resolve().parents[2]/'mujoco_arm/exports/anatomical_concept/elbow_guides'/('ground_frame_swept' if '--swept-root' in sys.argv else 'ground_frame_low_root' if '--low-root' in sys.argv else 'ground_frame_candidate')

if __name__=='__main__':
    audit.DEST=DEST;audit.run()
    p=DEST/'bone_motion_audit.json';a=json.loads(p.read_text(encoding='utf8'));a.update(scope='M7 known rigid skeleton, adapter-integrated frames, nominal key-mount hardware and sheaves; all included rigid pairs and bony references checked',not_included=['44 TPU/free distal instances, including two uncrimped old stops; no cap or muscle motion law','complete passive cable transmission and actuation','manufacturing tolerance, preload, material or loads']);p.write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf8')
