"""Interactive zero-pose arm/interface viewer. No actuation simulation."""
from pathlib import Path
import argparse,json,time,xml.etree.ElementTree as ET
import numpy as np
import mujoco

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'mujoco_arm/exports/skeleton_from_v2/V2-S1'
LAYOUT=ROOT/'mujoco_arm/exports/wire_psm_from_v2s1/PSM-W1/zero_scene.json'

def vec(x):return ' '.join(str(float(v)) for v in np.asarray(x).ravel())
def quat(matrix):
    q=np.empty(4);mujoco.mju_mat2Quat(q,np.asarray(matrix,dtype=float).ravel());return vec(q)

def build_model():
    root=ET.parse(BASE/'skeleton.xml').getroot();root.set('model','PSM-W1 zero-pose layout review')
    root.find('compiler').set('meshdir',str(BASE));root.find('option').set('gravity','0 0 0')
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='site' or (child.tag=='geom' and child.get('name')=='floor'):parent.remove(child)
    for g in root.findall('.//geom'):g.set('group','4' if g.get('name')=='thorax_geom' else '0')
    world=root.find('worldbody')
    def geom(name,kind,group,rgba,**attrs):
        ET.SubElement(world,'geom',{'name':name,'type':kind,'group':str(group),'rgba':vec(rgba),'contype':'0','conaffinity':'0','density':'0',**attrs})
    def line(name,a,b,rgba):
        if np.linalg.norm(np.asarray(b)-a)>1e-8:geom(name,'capsule',1,rgba,fromto=vec(np.r_[a,b]*.001),size='.0007')
    for z in json.loads(LAYOUT.read_text()):
        name=z['id'];a=np.array(z['A_mj_mm']);c=np.array(z['C_mj_mm']);length=np.linalg.norm(c-a)
        color=[.2,.72,.46,.3] if name.startswith('E') else [.2,.55,.9,.3] if name.startswith('G') else [.95,.6,.18,.3] if name.startswith('S') else [.8,.35,.62,.3]
        geom(name+'_belly_space','box',2,color,pos=vec((a+c)*.0005),size=vec([length*.0005,.017,.007]),quat=quat(z['belly_basis']))
        line(name+'_active',a,c,[.15,.35,.95,1]);p=np.array(z['passive_path_mj_mm'])
        for i,(start,end) in enumerate(zip(p[:-1],p[1:])):line(name+'_tendon_'+str(i),start,end,[.92,.46,.12,1])
        geom(name+'_sheave_envelope','cylinder',5,[.68,.7,.72,1],pos=vec(np.array(z['guide_center_mj_mm'])*.001),size='.0095 .003',quat=quat(z['guide_basis']))
        for role,color in [('A',[.1,.7,.2,1]),('C',[.1,.3,1,1]),('G',[1,.8,.1,1]),('B',[.9,.15,.2,1])]:geom(name+'_'+role,'sphere',3,color,pos=vec(np.array(z[role+'_mj_mm'])*.001),size='.0028')
    return mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true',help='Load and verify without opening a window');args=parser.parse_args()
    m=build_model();d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    if args.check:
        assert np.isfinite(d.geom_xpos).all();print({'geometries':m.ngeom,'bone_joint_coordinates':m.nq,'layout_modules':16,'view':'zero pose; no simulation'});return
    from mujoco import viewer as mjviewer
    pending=[]
    def key(code):pending.append(code)
    with mjviewer.launch_passive(m,d,key_callback=key) as view:
        def reset():
            view.cam.lookat[:]=[.07,0,0];view.cam.distance=1.3;view.cam.azimuth=130;view.cam.elevation=-10
        reset();view.opt.flags[mujoco.mjtVisFlag.mjVIS_JOINT]=False
        print('ZERO-POSE REVIEW, no actuation. Mouse: rotate / pan / zoom. Keys B=limb bones, W=wires, M=belly spaces, P=points, T=thorax, G=guides, R=reset camera. Close window to exit.',flush=True)
        groups={ord('B'):0,ord('W'):1,ord('M'):2,ord('P'):3,ord('T'):4,ord('G'):5}
        while view.is_running():
            with view.lock():
                while pending:
                    code=pending.pop(0)
                    if code in groups:view.opt.geomgroup[groups[code]]=1-view.opt.geomgroup[groups[code]]
                    elif code==ord('R'):reset()
                # Geometry overlays are zero-pose references. Deliberately do
                # not run dynamics or let joint sliders detach them from bones.
                d.qpos[:]=0;mujoco.mj_forward(m,d)
            view.sync();time.sleep(1/60)
if __name__=='__main__':main()
