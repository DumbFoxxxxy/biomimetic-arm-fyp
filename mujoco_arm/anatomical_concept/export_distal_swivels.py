"""Standalone M5 candidate rigid STEP parts, not validated print/force claims."""
from export_SC_cad import export_cad,ROOT
import sys
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15'/('swivel_terminal_distal_load' if '--distal-load' in sys.argv else 'swivel_terminal_candidate')
def run():
    if '--distal-load' in sys.argv:
        from export_swivel_stop_cad import run as native_run
        native_run();return
    for fid in ['FJ1486','FJ1477']:
        export_cad.OUT=BASE/fid;export_cad.run()
if __name__=='__main__':run()
