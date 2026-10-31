"""Actual MuJoCo mesh viewer; no muscle deformation/drive animation invented."""
from pathlib import Path
import json,sys,xml.etree.ElementTree as ET
import numpy as np,mujoco
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept/elbow_muscles'

def build(fid):
    folder=OUT/fid;data=json.loads((folder/'design.json').read_text(encoding='utf8'));root=ET.Element('mujoco',model=fid+'_unmounted_PSM_geometry');ET.SubElement(root,'compiler',meshdir='.',angle='degree');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='1000');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1');ET.SubElement(world,'light',pos='0 1 0',dir='0 -1 0',diffuse='.4 .4 .4')
    for name in data['components']:
        ET.SubElement(asset,'mesh',name=name,file=name+'.obj',scale='.001 .001 .001')
        color='.8 .12 .12 1' if name=='TPU_muscle' else '.2 .5 .7 1' if name.endswith('rigid_cap') else '.55 .6 .65 1'
        ET.SubElement(world,'geom',name=name,mesh=name,type='mesh',rgba=color,contype='0',conaffinity='0',group='0')
    ET.SubElement(asset,'mesh',name='source',file='../../assets/'+fid+'.stl',scale='.001 .001 .001');ET.SubElement(world,'geom',name='source',mesh='source',type='mesh',rgba='.75 .6 .6 0',group='1',contype='0',conaffinity='0')
    nodes=np.array(data['wire_polyline_display_mm'])/1000
    for i,(a,b) in enumerate(zip(nodes,nodes[1:])):ET.SubElement(world,'geom',name='wire_'+str(i),type='capsule',fromto=' '.join(map(str,np.r_[a,b])),size='.0003',rgba='.08 .08 .08 1',group='2',contype='0',conaffinity='0')
    ET.indent(root);path=folder/'module_view.xml';ET.ElementTree(root).write(path,encoding='utf8',xml_declaration=True)
    return path,data

def camera(cam,fid,data):
    nodes=np.array(data['wire_polyline_display_mm'])/1000;cam.lookat[:]=nodes.mean(0);cam.distance=data['body_span_mm']/1000*1.9;cam.azimuth=-90 if fid=='FJ1486' else 90;cam.elevation=-4

def run():
    if '--render' in sys.argv:
        from PIL import Image,ImageDraw
        panels=[]
        for fid in ['FJ1486','FJ1477']:
            path,data=build(fid);m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d);renderer=mujoco.Renderer(m,900,600);cam=mujoco.MjvCamera();camera(cam,fid,data);original=m.geom_rgba.copy();frames=[]
            for source in [True,False]:
                for g in range(m.ngeom):
                    name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g)
                    m.geom_rgba[g,3]=(1 if name=='source' else 0) if source else (0 if name=='source' else original[g,3])
                renderer.update_scene(d,camera=cam);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,600,65],fill='white');draw.text((12,12),fid+(' | source anatomy' if source else ' | generated manufacturing geometry'),fill='black');draw.text((12,35),'Not installed / deformed / actuated; nominal wire only',fill='black');frames.append(im)
            panel=Image.new('RGB',(1200,900),'white');panel.paste(frames[0],(0,0));panel.paste(frames[1],(600,0));panel.save(OUT/fid/'source_vs_module.png');panels.append(panel);renderer.close()
        combined=Image.new('RGB',(1200,1800),'white');combined.paste(panels[0],(0,0));combined.paste(panels[1],(0,900));combined.save(OUT/'shaped_muscle_comparison.png');print(OUT/'shaped_muscle_comparison.png');return
    fid=sys.argv[sys.argv.index('--id')+1] if '--id' in sys.argv else 'FJ1486';path,data=build(fid)
    m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    from mujoco import viewer
    import time
    with viewer.launch_passive(m,d) as v:
        camera(v.cam,fid,data)
        while v.is_running():v.sync();time.sleep(.03)

if __name__=='__main__':run()
