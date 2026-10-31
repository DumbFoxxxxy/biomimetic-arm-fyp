"""Static neutral M2 assembly only; no unverified soft motion is animated."""
from pathlib import Path
import json,sys,os,xml.etree.ElementTree as ET
import numpy as np,mujoco
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'

def build():
    data=json.loads((DEST/'assembly.json').read_text());root=ET.Element('mujoco',model='M2_neutral_proximal_mounts_only');ET.SubElement(root,'compiler',meshdir='.',angle='degree');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='1000');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1');ET.SubElement(world,'light',pos='0 1 .5',dir='0 -1 0',diffuse='.4 .4 .4')
    for r in data['instances']:
        key=r['key'];ET.SubElement(asset,'mesh',name=key,file=os.path.relpath(ROOT/r['file'],DEST).replace('\\','/'),scale='.001 .001 .001')
        color='.8 .12 .12 1' if r['kind']=='TPU' else '.2 .45 .7 1' if 'adapter' in key or 'rigid_cap' in key else '.88 .85 .78 1' if key in ['FJ3362','FJ3384','FJ3368','FJ3391','FJ3290'] else '.5 .6 .65 1'
        ET.SubElement(world,'geom',name=key,type='mesh',mesh=key,pos=' '.join(str(v/1000) for v in r['translation_display_mm']),rgba=color,contype='0',conaffinity='0')
    ET.indent(root);path=DEST/'neutral_view.xml';ET.ElementTree(root).write(path,encoding='utf8',xml_declaration=True);return path,data

def camera(cam,data):
    cam.lookat[:]=[-.05,0,.02];cam.distance=1.1;cam.azimuth=-90;cam.elevation=-4
    if '--detail' in sys.argv:
        fid=sys.argv[sys.argv.index('--id')+1] if '--id' in sys.argv else 'FJ1477';p=next(p for p in data['placements'] if p['id']==fid);cam.lookat[:]=np.array(p['proximal_pad_centre_mm'])/1000;cam.distance=.32;cam.azimuth=90 if fid=='FJ1477' else -90;cam.elevation=-10

def run():
    path,data=build();m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    if '--render' in sys.argv:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(m,1000,900);cam=mujoco.MjvCamera();camera(cam,data);frames=[]
        for azimuth in ([cam.azimuth-30,cam.azimuth+30] if '--detail' in sys.argv else [-90,90]):
            cam.azimuth=azimuth;renderer.update_scene(d,camera=cam);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,900,65],fill='white');draw.text((12,12),'M2 | actual neutral meshes | proximal mounts only',fill='black');draw.text((12,35),'Distal tendon transmission / TPU deformation / loads NOT validated',fill='black');frames.append(im)
        panel=Image.new('RGB',(1800,1000),'white');panel.paste(frames[0],(0,0));panel.paste(frames[1],(900,0));name='mount_detail_'+(sys.argv[sys.argv.index('--id')+1] if '--id' in sys.argv else 'FJ1477')+'.png' if '--detail' in sys.argv else 'mounted_neutral.png';panel.save(DEST/name);renderer.close();print(DEST/name);return
    from mujoco import viewer
    import time
    with viewer.launch_passive(m,d) as v:
        camera(v.cam,data)
        while v.is_running():v.sync();time.sleep(.03)

if __name__=='__main__':run()
