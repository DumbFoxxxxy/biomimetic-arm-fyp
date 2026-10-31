"""Render the actual exported mechanism using the system Pillow environment."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json
import numpy as np
import mujoco
from PIL import Image,ImageDraw,ImageFont

EXPORT=Path(__file__).resolve().parent.parent/'exports/blade_motion_design'

def main():
    model=_load_project_model(str(EXPORT/'gimbal_arm.xml'));states={s['action']:s['qpos'] for s in json.loads((EXPORT/'gimbal_render_states.json').read_text())};renderer=mujoco.Renderer(model,height=480,width=640);panels=[]
    labels=['零位｜肩胛—同心三轴—肘部','联合抬臂 170°｜目标姿态','GH 零位｜内侧轴承、外侧开口','GH 外展 120°｜上臂通过开口']
    for label,action,detail in [(labels[0],'ZERO',False),(labels[1],'RAISE_170',False),(labels[2],'ZERO',True),(labels[3],'GH_ABD',True)]:
        d=mujoco.MjData(model);d.qpos[:]=states[action];mujoco.mj_forward(model,d);camera=mujoco.MjvCamera();camera.azimuth=-100;camera.elevation=-12
        if detail:
            bid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'gh_flex');camera.lookat[:]=d.xpos[bid];camera.distance=.35;camera.azimuth=-125;camera.elevation=-20
        else:camera.lookat[:]=[.13,-.02,.15];camera.distance=1.55
        option=mujoco.MjvOption();option.sitegroup[:]=0;option.geomgroup[3]=0;renderer.update_scene(d,camera=camera,scene_option=option);panels.append(Image.fromarray(renderer.render()))
    renderer.close();page=Image.new('RGB',(1280,1120),'#f3f5f8');draw=ImageDraw.Draw(page);font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',23);small=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',19)
    for i,panel in enumerate(panels):
        x=(i%2)*640;y=(i//2)*535;draw.text((x+18,y+14),labels[i],font=font,fill='#122238');page.paste(panel,(x,y+50))
    draw.text((18,1084),'运动设计候选；未制造。橙色线为钢丝，图中姿态由模型渲染。',font=small,fill='#334155');page.save(EXPORT/'gimbal_candidate_views.png');print(EXPORT/'gimbal_candidate_views.png')
if __name__=='__main__':main()
