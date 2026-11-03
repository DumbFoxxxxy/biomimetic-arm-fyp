"""Actual MuJoCo static view of the eight controlled axis/channel specimens."""
import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'mujoco_arm/tools')]
from project_runtime import load_mujoco_model
import build_shoulder_axis_coupons as C
import mesh_geom


def run(interactive=False):
    import mujoco
    source=C.DEST/'manifest.json';report=json.loads(source.read_text(encoding='utf8'))
    root=ET.Element('mujoco',model='STATIC_axis_coupons_NOT_calibrated')
    ET.SubElement(root,'compiler',angle='radian')
    visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='900')
    assets=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody')
    ET.SubElement(world,'light',pos='0 -.2 1',dir='0 0 -1')
    low=[];high=[]
    for index,coupon in enumerate(report['coupons']):
        offset=np.array([(index%4)*55.,(index//4)*55.,0.])
        for key,row in coupon['sources'].items():
            path=ROOT/row['file']
            if C.A.hashed(path)!=row['sha256']:raise ValueError('Changed coupon part '+str(path))
            name='coupon_'+str(index)+'_'+key
            ET.SubElement(assets,'mesh',name=name,file=str(path),scale='.001 .001 .001')
            pos=(offset+np.array(row['translation_source_mm']))*.001
            colour=('.85 .12 .10 1' if index<4 else '.1 .55 .85 1') if key=='TPU' else '.7 .7 .65 1'
            ET.SubElement(world,'geom',name=name,type='mesh',mesh=name,pos=' '.join(map(str,pos)),
                          rgba=colour,contype='0',conaffinity='0')
            v,_=mesh_geom.read_stl(str(path));low.append(v.min(0)+offset);high.append(v.max(0)+offset)
    scene=C.DEST/'comparison_scene.xml';ET.ElementTree(root).write(scene,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(scene);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    look=(np.min(low,axis=0)+np.max(high,axis=0))*.0005
    def camera(cam):
        cam.lookat[:]=look;cam.distance=.39;cam.azimuth=115;cam.elevation=30
    if interactive:
        from mujoco import viewer
        import time
        with viewer.launch_passive(model,data) as view:
            camera(view.cam)
            while view.is_running():view.sync();time.sleep(.03)
    else:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(model,height=900,width=1200);cam=mujoco.MjvCamera();camera(cam)
        renderer.update_scene(data,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        image=Image.fromarray(renderer.render());draw=ImageDraw.Draw(image)
        draw.rectangle((0,0,1200,80),fill='white')
        draw.text((12,10),'STATIC AXIS / CHANNEL COUPONS | RED: 15/10/30 | BLUE: 10/8/10',fill='black')
        draw.text((12,33),'Each set: Y bare, Y tube, Z bare, Z tube. Same XYZ print axes; only two opposed grip plates.',fill='black')
        draw.text((12,55),'No measured force or stroke. Tube variants show nominal washer, uncrimped stop and wire route.',fill='black')
        image.save(C.DEST/'comparison.png');renderer.close()
    (C.DEST/'view_check.json').write_text(json.dumps(dict(manifest_sha256=C.A.hashed(source),
        script_sha256=C.A.hashed(Path(__file__).resolve()),geometry_instances=model.ngeom,DOF=model.nq,
        static_only=True,stage_pass=False),indent=2)+'\n',encoding='utf8',newline='\n')
    print('Coupon scene',model.ngeom,'objects, DOF',model.nq)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--view',action='store_true')
    run(parser.parse_args().view)
