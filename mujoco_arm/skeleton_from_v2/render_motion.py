"""Render prescribed accepted paths with MuJoCo; no dynamic/force claim."""
from pathlib import Path
import json,numpy as np,mujoco
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parents[1]/'exports/skeleton_from_v2'
def run():
    m=mujoco.MjModel.from_xml_path(str(OUT/'skeleton_review.xml'));d=mujoco.MjData(m);renderer=mujoco.Renderer(m,480,640)
    camera=mujoco.MjvCamera();camera.lookat[:]=[.08,0,.23];camera.distance=2.2;camera.azimuth=130;camera.elevation=-10
    paths=json.loads((OUT/'retained_shape_paths.json').read_text());frames=[]
    for name,seq in paths.items():
        indexes=np.unique(np.linspace(0,len(seq)-1,18).astype(int));forward=[seq[i] for i in indexes]
        for q in forward+forward[-2::-1]:
            d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d);renderer.update_scene(d,camera=camera)
            im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,640,42],fill='white');draw.text((12,8),'V2-S1 | '+name+' | prescribed bone motion; no actuation',fill='black');frames.append(im)
        if forward:frames[-len(forward)].save(OUT/(name+'_endpoint.png'))
    renderer.close();frames[0].save(OUT/'skeleton_motion_review.gif',save_all=True,append_images=frames[1:],duration=70,loop=0)
    print('frames',len(frames),flush=True)
if __name__=='__main__':run()
