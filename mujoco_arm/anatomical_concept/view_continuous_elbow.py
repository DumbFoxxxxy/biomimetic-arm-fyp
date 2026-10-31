"""Actual continuous-channel elbow, neutral geometry only. No TPU animation."""
from pathlib import Path
import json, sys, xml.etree.ElementTree as ET
import numpy as np
import mujoco

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
DEST = BASE/'continuous_channel_candidate'


def run():
    tree = ET.parse(BASE/'neutral_view.xml')
    for fid in ['FJ1486', 'FJ1477']:
        mesh = tree.find(f"./asset/mesh[@name='{fid}_TPU_muscle']")
        assert mesh is not None
        mesh.set('file',str(DEST/fid/'TPU_muscle.obj'))
    tree.getroot().set('model','continuous_channel_elbow_neutral_UNVALIDATED')
    tree.write(DEST/'neutral_view.xml', encoding='utf8', xml_declaration=True)
    m=mujoco.MjModel.from_xml_path(str(DEST/'neutral_view.xml'));d=mujoco.MjData(m)
    mujoco.mj_forward(m,d)
    print('neutral scene compiled; geoms',m.ngeom,'DOF',m.nq,flush=True)
    if '--check' in sys.argv: return
    def camera(cam):
        cam.lookat[:]=[-.07,0,-.03];cam.distance=.61;cam.azimuth=125;cam.elevation=-7
    if '--render' in sys.argv:
        from PIL import Image, ImageDraw
        m.vis.global_.offwidth=1400;m.vis.global_.offheight=1000
        renderer=mujoco.Renderer(m,height=1000,width=1400);cam=mujoco.MjvCamera();camera(cam)
        renderer.update_scene(d,camera=cam)
        renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im)
        draw.rectangle((0,0,1400,65),fill='white')
        draw.text((15,12),'ACTUAL ANATOMICAL ELBOW + SHAPED TPU LATTICES + CONTINUOUS WIRE CHANNELS',fill='black')
        draw.text((15,36),'Neutral manufacturing candidate | no verified TPU deformation, drive or physical performance',fill='black')
        im.save(DEST/'elbow_neutral.png');renderer.close();return
    import time
    from mujoco import viewer
    with viewer.launch_passive(m,d) as v:
        camera(v.cam)
        while v.is_running(): v.sync();time.sleep(.03)


if __name__=='__main__':run()
