"""Eleven parts of the connected paired-cap support, nominal neutral pose."""
from export_SC_cad import export_cad,ROOT
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support'
if __name__=='__main__':export_cad.run()
