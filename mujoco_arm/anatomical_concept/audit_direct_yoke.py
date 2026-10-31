"""All-pair rigid skeleton/interface paths for isolated direct-muscle candidate."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
import sys,json
from project_runtime import WorkspacePath as Path
sys.argv.append('--muscle-mounts')
import audit_joint_model as audit
DEST=Path(__file__).resolve().parents[2]/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces'/('mechanical_yoke_shifted' if '--shifted' in sys.argv else 'mechanical_yoke_zero')
if '--access-relief' in sys.argv:DEST=DEST.parent/'yoke_access_relief_study'
if '--extensor-reanchor' in sys.argv:DEST=DEST.parent/'mechanical_yoke_reanchored_pair'
if __name__=='__main__':
    audit.DEST=DEST;audit.run();p=DEST/'bone_motion_audit.json';a=json.loads(p.read_text(encoding='utf8'));a.update(scope='M9 direct-muscle mechanical end constraints:103 nominal rigid attached pieces, all pairs/bone references; no TPU strain or actuation proof',not_included=['2TPU flexible meshes and2uncrimped cable-stop nominal poses','muscle deformation/wire paths during motion','strength,crimp,manufacturing/material calibration ordrive']);p.write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf8')
