"""Reuse generic faceted-BRep STEP exporter, without inheriting old arm geometry."""
from pathlib import Path
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
