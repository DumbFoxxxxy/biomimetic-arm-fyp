"""Direct prototype six new/revised rigid parts, faceted editable STEP."""
from export_SC_cad import export_cad,ROOT
import sys
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces'/('mechanical_yoke_shifted' if '--shifted' in sys.argv else 'mechanical_yoke_zero')
if '--access-relief' in sys.argv:export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/yoke_access_relief_study'
if '--extensor-reanchor' in sys.argv:export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair'
if __name__=='__main__':export_cad.run()
