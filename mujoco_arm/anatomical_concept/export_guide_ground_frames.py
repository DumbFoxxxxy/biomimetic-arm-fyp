"""Standalone M7 ten new/revised rigid manufacturing parts."""
from export_SC_cad import export_cad, ROOT
import sys
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_guides'/('ground_frame_swept' if '--swept-root' in sys.argv else 'ground_frame_low_root' if '--low-root' in sys.argv else 'ground_frame_candidate')
if __name__=='__main__':export_cad.run()
