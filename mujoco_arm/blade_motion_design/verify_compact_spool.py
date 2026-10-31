"""Search wrap-side choices for the compact GR fairlead geometry."""
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from design_axis_spool import BASE,build_stable_spool
from shoulder_model import EXPORT,v8

def evaluate(s1,s2,n=161):
    p={**BASE,'radius_mm':25.,'entry_radius_mm':3.,'entry_width_mm':2.,'entry_guide_mm':[0.,19.,21.],'external_anchor_mm':[0.,52.,13.],'sector_side_angle_E01_deg':s1,'sector_side_angle_E02_deg':s2}
    root=build_stable_spool(p).getroot()
    for rid in ['E01','E02']:
        e=root.find(f'.//geom[@name="entry_{rid}"]');b=root.find(f'.//site[@name="bridge_{rid}"]');x=np.fromstring(e.get('pos'),sep=' ');y=np.fromstring(b.get('pos'),sep=' ');b.set('pos',v8.vec(x+(y-x)*8/25))
    for node in root.find('tendon'):
        for item in list(node)[:2]:node.remove(item)
    m=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);J=[];L=[]
    for theta in np.linspace(0,160,n):
        d.qpos[0]=np.radians(theta);mujoco.mj_forward(m,d);J.append(d.ten_J.copy());L.append(d.ten_length.copy())
    J=np.array(J);L=np.array(L);okay=bool(np.all(J[:,0]*J[:,1]<-1e-8) and np.min(np.abs(J))>.0249 and np.max(abs(np.diff(L,axis=0)))<.025*np.radians(160/(n-1))*1.1)
    return {'pass':okay,'sides_deg':[s1,s2],'samples':n,'min_arm_mm':float(np.min(np.abs(J))*1000),'max_length_step_mm':float(abs(np.diff(L,axis=0)).max()*1000),'parameters':p}

def main():
    original=evaluate(330.,0.,2561);print('ORIGINAL',original,flush=True);best=original
    if not best['pass']:
        for s1 in np.arange(0,360,15):
            for s2 in np.arange(0,360,15):
                row=evaluate(float(s1),float(s2),161)
                if row['pass']:
                    row=evaluate(float(s1),float(s2),2561)
                    if row['pass']:best=row;break
            if best['pass']:break
    (EXPORT/'compact_spool_verification.json').write_text(json.dumps(best,indent=2));print('BEST',best,flush=True)
if __name__=='__main__':main()
