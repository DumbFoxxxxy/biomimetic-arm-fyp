"""Independent faceted engineering CAD of two real receiver mechanisms."""
from export_SC_cad import export_cad,ROOT
import sys
if __name__=='__main__':
    if '--mount-only' in sys.argv:
        export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/forearm_receivers'/('rephased_anchors' if '--rephase' in sys.argv else '');export_cad.run()
    else:
        for fid in ['FJ1486','FJ1477']:
            export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/forearm_receivers'/fid;export_cad.run()
