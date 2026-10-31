from pathlib import Path
import xml.etree.ElementTree as ET
import mujoco
from PIL import Image,ImageDraw,ImageFont

p=Path(__file__).resolve().parents[1]/'exports/tpu_muscle_design'
root=ET.Element('mujoco');asset=ET.SubElement(root,'asset');ET.SubElement(asset,'mesh',name='coupon',file=str(p/'psm_reference_coupon_mm.stl'),scale='.001 .001 .001')
world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='.05 -.1 .2');ET.SubElement(world,'geom',type='mesh',mesh='coupon',rgba='.1 .7 .6 1');ET.SubElement(world,'geom',type='capsule',fromto='-.025 0 0 .100 0 0',size='.0005',rgba='1 .4 .05 1')
m=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);r=mujoco.Renderer(m,height=480,width=640);c=mujoco.MjvCamera();c.lookat[:]=[.0375,0,0];c.distance=.19;c.azimuth=-110;c.elevation=-24;r.update_scene(d,camera=c);im=Image.fromarray(r.render());r.close()
page=Image.new('RGB',(640,565),'#f4f6f8');page.paste(im,(0,50));g=ImageDraw.Draw(page);font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',20);g.text((12,12),'TPU 晶格内部结构样件｜尚未验证收缩',font=font,fill='#123');g.text((12,535),'单胞段 75 mm；绿色 TPU；橙色为内部钢丝',font=font,fill='#123');page.save(p/'coupon_preview.png')
