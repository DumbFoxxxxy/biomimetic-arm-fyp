"""Export only render states; this does not freeze or approve a design."""
import json,xml.etree.ElementTree as ET
import mujoco
from gimbal_arm_model import build_gimbal
from arm_targets import target_actions
from audit_gimbal_drive import apply_target
from shoulder_model import EXPORT

if __name__=='__main__':
    tree,_=build_gimbal();tree.write(EXPORT/'gimbal_arm.xml',encoding='utf-8',xml_declaration=True);m=mujoco.MjModel.from_xml_string(ET.tostring(tree.getroot(),encoding='unicode'));d=mujoco.MjData(m);actions=target_actions();out=[]
    for name,z in [('ZERO',actions['GH_FLEX'][0]),('RAISE_170',actions['RAISE_170'][-1]),('GH_ABD',actions['GH_ABD'][-1])]:
        apply_target(m,d,z);out.append({'action':name,'qpos':d.qpos.tolist()})
    (EXPORT/'gimbal_render_states.json').write_text(json.dumps(out,indent=2))
