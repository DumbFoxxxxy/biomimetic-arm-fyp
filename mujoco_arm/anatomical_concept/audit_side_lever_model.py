"""Reuse all-pair prescribed motion probe on isolated M4 sources."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15'
    folder=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate'/variant
    if '--neutral' in sys.argv:
        import check_elbow_muscle_mounts as check
        check.DEST=folder;check.run();return
    if '--exports' in sys.argv:
        import check_mount_exports as check
        check.DEST=folder;check.run();return
    sys.argv.append('--muscle-mounts')
    import audit_joint_model as audit
    audit.DEST=folder;audit.run()
if __name__=='__main__':run()
