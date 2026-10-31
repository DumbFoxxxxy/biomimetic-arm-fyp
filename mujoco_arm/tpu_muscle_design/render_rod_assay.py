from pathlib import Path
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from PIL import Image,ImageDraw,ImageFont
OUT=Path(__file__).resolve().parents[1]/'exports/tpu_muscle_design'
table=np.load(OUT/'coupon_rod_equilibria.npz');branches=table['branches'];panels=[]
for index in [0,10,12]:
    points=table['states'][index]*.001;span=float(points[:,0].max());root=ET.Element('mujoco');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='.04 -.15 .2')
    for line in branches:
        for a,b in zip(line[:-1],line[1:]):ET.SubElement(world,'geom',type='capsule',fromto=' '.join(map(str,np.r_[points[a],points[b]])),size='.001',rgba='.12 .7 .6 1')
    ET.SubElement(world,'geom',type='capsule',fromto=f'-.02 0 0 {span+.02} 0 0',size='.0005',rgba='1 .45 .1 1')
    for x in [-.005,span+.005]:ET.SubElement(world,'geom',type='box',pos=f'{x} 0 0',size='.005 .006 .016',rgba='.12 .5 .45 1')
    model=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));data=mujoco.MjData(model);mujoco.mj_forward(model,data);r=mujoco.Renderer(model,height=400,width=640);c=mujoco.MjvCamera();c.lookat[:]=[span/2,0,0];c.distance=.18;c.azimuth=-110;c.elevation=-24;r.update_scene(data,camera=c);panels.append(Image.fromarray(r.render()));r.close()
page=Image.new('RGB',(1920,505),'#f4f6f8');draw=ImageDraw.Draw(page);font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',22)
for i,(im,label) in enumerate(zip(panels,['零位｜单胞段 75 mm','压缩 25%｜杆模型平衡构形','压缩 30%｜已检出接触，当前无接触模型失效'])):
    page.paste(im,(i*640,50));draw.text((i*640+12,12),label,font=font,fill='#123')
draw.text((15,466),'离散杆计算示意：非实测；未包含外层、真实 TPU 本构、导向环刚度及接触反力。',font=font,fill='#123');page.save(OUT/'coupon_rod_assay.png')
