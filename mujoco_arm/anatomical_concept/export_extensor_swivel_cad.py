"""Six local swivel candidates only; no bone-grounded support claimed."""
from export_SC_cad import export_cad,ROOT
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/cap_swivel_study'
if __name__=='__main__':export_cad.run()
