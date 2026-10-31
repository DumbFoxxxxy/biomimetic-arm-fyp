"""Render exported XML with the system Python (MuJoCo and Pillow installed)."""
from pathlib import Path
import json
import numpy as np
import mujoco
from PIL import Image,ImageDraw,ImageFont

EXPORT=Path(__file__).resolve().parent.parent/'exports/blade_motion_design'
def render():
    model=mujoco.MjModel.from_xml_path(str(EXPORT/'shoulder_no_contact.xml'))
    poses=json.loads((EXPORT/'render_states.json').read_text())
    renderer=mujoco.Renderer(model,height=480,width=640)
    camera=mujoco.MjvCamera();camera.lookat[:]=[.13,-.035,.19];camera.distance=.64;camera.azimuth=-100;camera.elevation=-8
    panels=[]
    for pose in poses:
        data=mujoco.MjData(model);data.qpos[:]=pose['qpos'];mujoco.mj_forward(model,data)
        renderer.update_scene(data,camera=camera);panels.append(Image.fromarray(renderer.render()).convert('RGB'))
    renderer.close()
    page=Image.new('RGB',(1280,570),'#f4f5f7');draw=ImageDraw.Draw(page)
    font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',22)
    for index,panel in enumerate(panels):page.paste(panel,(640*index,55))
    draw.text((25,15),'拟合零位：八根钢丝＋S07局部通道',font=font,fill='#16243b')
    draw.text((665,15),'肩胛上回旋50°：保留SC/AC关节中心',font=font,fill='#16243b')
    draw.text((25,538),'几何与无载仿真候选；并非已制造样机。钢丝固定端需机架支撑。',font=font,fill='#334155')
    page.save(EXPORT/'shoulder_candidate_views.png')
    print(EXPORT/'shoulder_candidate_views.png')
if __name__=='__main__':render()
