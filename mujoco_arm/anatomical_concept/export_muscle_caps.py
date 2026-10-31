"""Editable faceted STEP for rigid terminal caps; TPU remains a print mesh."""
from export_SC_cad import export_cad,ROOT
import json
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles'

def run():
    for fid in ['FJ1486','FJ1477']:
        folder=BASE/fid;data=json.loads((folder/'cad_input.json').read_text(encoding='utf8'));dest=folder/'hard_caps';dest.mkdir(exist_ok=True)
        parts=[p for p in data['parts'] if p['name'].endswith('rigid_cap')]
        (dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':parts}),encoding='utf8')
        export_cad.OUT=dest;export_cad.run()

if __name__=='__main__':run()
