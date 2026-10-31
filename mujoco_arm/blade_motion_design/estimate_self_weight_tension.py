"""Report the unbounded quasistatic minimax force, without choosing hardware."""
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from scipy.optimize import linprog
from gimbal_arm_model import build_gimbal
from arm_targets import target_actions
from audit_gimbal_drive import apply_target,native_jacobian
from shoulder_model import EXPORT
from verification_snapshot import fingerprint

if __name__=='__main__':
    tree,_=build_gimbal();root=tree.getroot();root.find('option').set('gravity','0 0 -9.81');m=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);rows=[]
    for action,seq in target_actions().items():
        peak=0.;worst=None
        for i,z in enumerate([a*(1-f)+b*f for a,b in zip(seq[:-1],seq[1:]) for f in np.linspace(0,1,5,endpoint=False)]+[seq[-1]]):
            apply_target(m,d,z);d.qvel[:]=0.;d.qacc[:]=0.;mujoco.mj_inverse(m,d);J=native_jacobian(m,d)
            lp=linprog(np.r_[np.zeros(m.ntendon),1.],A_ub=np.column_stack([np.eye(m.ntendon),-np.ones(m.ntendon)]),b_ub=np.zeros(m.ntendon),A_eq=np.column_stack([J.T,np.zeros(m.nv)]),b_eq=-d.qfrc_inverse,bounds=[(.5,None)]*(m.ntendon+1),method='highs')
            if not lp.success:raise RuntimeError('No positive self-weight equilibrium')
            if lp.x[-1]>peak:
                peak=float(lp.x[-1]);worst={'pose_index':i,'largest_force_line':mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_TENDON,int(np.argmax(lp.x[:-1])))}
        row={'action':action,'minimum_required_peak_line_force_N':peak,**worst};rows.append(row);print(row,flush=True)
    report={'scope':'Estimated uniform 1100 kg/m3 solid parts; static gravity only, no payload, no friction; not a motor or cable rating. The unbounded optimization estimates required peak force instead of increasing any actuator limit.','actions':rows,'maximum_N':max(r['minimum_required_peak_line_force_N'] for r in rows),'input_fingerprint':fingerprint(root)}
    (EXPORT/'self_weight_tension_estimate.json').write_text(json.dumps(report,indent=2))
