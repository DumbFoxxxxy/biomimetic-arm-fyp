"""Individual M4 candidate rigid STEP exports, no drive validation."""
from export_SC_cad import export_cad,ROOT
import sys
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15'
    export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate'/variant
    export_cad.run()
if __name__=='__main__':run()
