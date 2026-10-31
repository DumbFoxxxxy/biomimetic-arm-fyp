"""MuJoCo exploded assembly visualization of the checked E1 installation order."""
from pathlib import Path
import json,xml.etree.ElementTree as ET
import numpy as np,mujoco
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'

def run():
    data=json.loads((OUT/'assembly.json').read_text(encoding='utf8'));checks=json.loads((OUT/'assembly_and_retention_checks.json').read_text(encoding='utf8'))
    if not checks['sequential_installation_geometry_pass']:raise RuntimeError('No checked assembly sequence')
    names=list(data['components']);axis=np.array(data['axis_source']);direction=np.array(checks['chosen_cartridge_path']['direction']);E=np.array(data['centre_source_mm'])-[-140,-70,1100]
    root=ET.Element('mujoco',model='E1_geometry_only_assembly');ET.SubElement(root,'compiler',meshdir=str(OUT),angle='degree');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='900',offheight='800');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1')
    for name in names:
        ET.SubElement(asset,'mesh',name=name,file=name+'.stl',scale='.001 .001 .001')
        body=ET.SubElement(world,'body',name=name,mocap='true');ET.SubElement(body,'geom',name=name,type='mesh',mesh=name,contype='0',conaffinity='0',rgba='.88 .85 .78 1' if name in ['humerus_GH_elbow_fork','forearm_elbow_hub'] else '.3 .55 .65 1')
    ET.indent(root);path=OUT/'assembly_view.xml';ET.ElementTree(root).write(path,encoding='utf8',xml_declaration=True)
    m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);bodies={n:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,n) for n in names};gids={n:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,n) for n in names}
    renderer=mujoco.Renderer(m,800,900);cam=mujoco.MjvCamera();cam.lookat[:]=E/1000;cam.distance=.38;cam.azimuth=-75;cam.elevation=-10
    stages=[('Seat forearm + bush',['forearm_elbow_hub','elbow_bearing_bush'],direction*.045,28),('Insert left sleeve',['elbow_spacer_left'],-axis*.085,20),('Insert right sleeve',['elbow_spacer_right'],axis*.085,20),('Insert pin + head-side washer',['elbow_shoulder_pin','elbow_washer_left'],-axis*(checks['pin_insertion']['distance_mm']/1000),28),('Position far washer',['elbow_washer_right'],axis*.025,16),('Position nut envelope',['elbow_retaining_nut'],axis*.025,16)]
    frames=[];installed={'humerus_GH_elbow_fork'}
    for title,moving,offset,count in stages:
        for i in range(count):
            f=i/(count-1)
            for name in names:
                visible=name in installed or name in moving;m.geom_rgba[gids[name],3]=int(visible)
                d.mocap_pos[m.body_mocapid[bodies[name]]]=offset*(1-f) if name in moving else 0
            mujoco.mj_forward(m,d);renderer.update_scene(d,camera=cam);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,900,65],fill='white');draw.text((15,12),title+' | E1 nominal assembly geometry',fill='black');draw.text((15,35),'Threads / preload / printed fit / strength NOT validated',fill='black');frames.append(im)
        installed.update(moving)
    frames[0].save(OUT/'elbow_assembly.gif',save_all=True,append_images=frames[1:],duration=90,loop=0);frames[-1].save(OUT/'elbow_assembled.png');renderer.close();print('assembly frames',len(frames))

if __name__=='__main__':run()
