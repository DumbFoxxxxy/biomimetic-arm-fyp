"""M6 cassette faceted engineering STEP; no mates or feature-history claim."""
from export_SC_cad import export_cad, ROOT
export_cad.OUT=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_guides/template'
if __name__=='__main__':export_cad.run()
