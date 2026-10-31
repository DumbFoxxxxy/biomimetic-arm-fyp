"""Render the coordinate proposal; boxes are space reservations, not TPU FEM."""
from pathlib import Path
import json,numpy as np,mujoco
from PIL import Image,ImageDraw,ImageFont
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'mujoco_arm/exports/wire_psm_from_v2s1';DEST=OUT/'PSM-W1'
def run():
    m=mujoco.MjModel.from_xml_path(str(ROOT/'mujoco_arm/exports/skeleton_from_v2/V2-S1/skeleton.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);m.vis.global_.offwidth=640;m.vis.global_.offheight=640
    renderer=mujoco.Renderer(m,640,640);modules=json.loads((DEST/'zero_scene.json').read_text());panels=[]
    def add(scene,kind,size,pos,mat,color):
        if scene.ngeom>=scene.maxgeom:raise RuntimeError('Rendering capacity exceeded')
        g=scene.geoms[scene.ngeom];mujoco.mjv_initGeom(g,kind,np.array(size),np.array(pos),np.array(mat).ravel(),np.array(color));scene.ngeom+=1
    def line(scene,a,b,color,radius=.0008):
        a=np.array(a)*.001;b=np.array(b)*.001;delta=b-a;length=np.linalg.norm(delta)
        if length<1e-10:return
        axis=delta/length;helper=[1,0,0] if abs(axis[0])<.8 else [0,1,0];x=np.cross(helper,axis);x/=np.linalg.norm(x);y=np.cross(axis,x);add(scene,mujoco.mjtGeom.mjGEOM_CAPSULE,[radius,length/2,0],(a+b)/2,np.column_stack([x,y,axis]),color)
    for view in ['all','E01']:
        camera=mujoco.MjvCamera();camera.lookat[:]=[.07,0,0];camera.distance=1.25;camera.azimuth=130;camera.elevation=-10
        if view=='E01':camera.lookat[:]=[.18,0,.02];camera.distance=.75;camera.azimuth=15
        if view=='E01':
            for name in ['thorax_geom','clavicle_visual','scapula_mesh']:
                gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,name);m.geom_group[gid]=3
        option=mujoco.MjvOption();option.geomgroup[3]=0
        renderer.update_scene(d,camera=camera,scene_option=option);scene=renderer.scene
        for module in modules:
            if view!='all' and module['id']!=view:continue
            A=np.array(module['A_mj_mm']);C=np.array(module['C_mj_mm']);G=np.array(module['G_mj_mm']);B=np.array(module['B_mj_mm']);basis=np.array(module['belly_basis']);s=np.linalg.norm(C-A)
            color=[.2,.72,.46,.35] if module['id'].startswith('E') else [.2,.55,.9,.35] if module['id'].startswith('G') else [.95,.6,.18,.3] if module['id'].startswith('S') else [.8,.35,.62,.3]
            add(scene,mujoco.mjtGeom.mjGEOM_BOX,[s*.0005,.017,.007],(A+C)*.0005,basis,color);line(scene,A,C,[.15,.35,.95,1]);p=np.array(module['passive_path_mj_mm'])
            for a,b in zip(p[:-1],p[1:]):line(scene,a,b,[.92,.46,.12,1])
            add(scene,mujoco.mjtGeom.mjGEOM_CYLINDER,[.0095,.003,0],np.array(module['guide_center_mj_mm'])*.001,module['guide_basis'],[.68,.7,.72,1])
            for point,c in [(A,[.1,.7,.2,1]),(C,[.1,.3,1,1]),(G,[1,.8,.1,1]),(B,[.9,.15,.2,1])]:add(scene,mujoco.mjtGeom.mjGEOM_SPHERE,[.0028,0,0],point*.001,np.eye(3),c)
        panels.append(Image.fromarray(renderer.render()))
    renderer.close();image=Image.new('RGB',(1280,740),'white');image.paste(panels[0],(0,55));image.paste(panels[1],(640,55));draw=ImageDraw.Draw(image)
    font=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',19);small=ImageFont.truetype('C:/Windows/Fonts/msyh.ttc',16)
    draw.text((20,14),'PSM-W1：全臂接口布局候选',fill='black',font=font);draw.text((660,14),'E01：肘部接口（仅显示肱骨和尺骨）',fill='black',font=font)
    draw.rectangle([0,650,1280,740],fill='white');draw.text((20,658),'绿色点 A：固定端　蓝色点 C：移动端　黄色点 G：进线切点　红色点 B：骨骼止点',fill='black',font=small)
    draw.text((20,687),'蓝线为肌腹内部主动钢丝，橙线为被动肌腱；半透明盒是安装空间情景，不是晶格变形或驱动验证。',fill='black',font=small)
    image.save(DEST/'layout_review.png');print('layout image saved',flush=True)
if __name__=='__main__':run()
