"""Reuse generic faceted-BRep STEP exporter, without inheriting old arm geometry."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'mujoco_arm/anatomical_arm'))
import export_cad
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection'
if '--AC' in sys.argv:export_cad.OUT=export_cad.OUT/('AC_compact' if '--compact' in sys.argv else 'AC_connection')
if '--GH' in sys.argv:export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact'/('GH_closed' if '--closed' in sys.argv else 'GH_preserve_head' if '--preserve-head' in sys.argv else 'GH_connection')
if '--elbow' in sys.argv:export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
if '--muscle-mounts' in sys.argv:export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'
if __name__=='__main__':export_cad.run()
