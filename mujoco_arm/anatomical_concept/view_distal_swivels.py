"""Actual M5 caps and local articulation preview, not TPU actuation."""
from pathlib import Path
import os,json,sys,hashlib,xml.etree.ElementTree as ET
import numpy as np,mujoco
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15'/('swivel_terminal_distal_load' if '--distal-load' in sys.argv else 'swivel_terminal_candidate')
def run():
    assembly=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));inputs={str((ROOT/r['file']).resolve()):hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest() for r in assembly['instances']}
    for name in ['assembly.json','FJ1486/design.json','FJ1477/design.json']:inputs[str((DEST/name).resolve())]=hashlib.sha256((DEST/name).read_bytes()).hexdigest()
    root=ET.Element('mujoco',model='M5_neutral_real_terminal_geometry');ET.SubElement(root,'compiler',meshdir='.',angle='degree');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1000',offheight='850');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='1 -1 .5',dir='-1 1 -1');ET.SubElement(world,'light',pos='-1 1 .5',dir='1 -1 -1',diffuse='.5 .5 .5')
    for r in assembly['instances']:
        key=r['key'];ET.SubElement(asset,'mesh',name=key,file=os.path.relpath(ROOT/r['file'],DEST).replace('\\','/'),scale='.001 .001 .001')
        color='.9 .6 .1 1' if key.endswith('ball_crimp_blank') else '.15 .5 .8 1' if r.get('m5_part') in ['distal_socket_base','socket_retainer'] else '.8 .12 .1 1' if r['kind']=='TPU' else '.85 .82 .75 1' if key in ['FJ3362','FJ3384','FJ3368','FJ3391','FJ3290'] else '.5 .6 .65 1'
        if '--cutaway' in sys.argv and key.endswith('socket_retainer'):color='.15 .5 .8 .28'
        parent=ET.SubElement(world,'body',name=key+'_motion',mocap='true') if key.endswith('ball_crimp_blank') else world
        ET.SubElement(parent,'geom',name=key,type='mesh',mesh=key,pos=' '.join(str(c/1000) for c in r['translation_display_mm']),rgba=color,contype='0',conaffinity='0')
    ET.indent(root);path=DEST/('neutral_cutaway_view.xml' if '--cutaway' in sys.argv else 'neutral_view.xml');ET.ElementTree(root).write(path,encoding='utf8',xml_declaration=True);m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    # Rendering-only isolation avoids a forearm lever blocking the small socket.
    # All131 meshes remain in the model and neutral collision/export checks.
    for r in assembly['instances']:
        gid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,r['key']);m.geom_group[gid]=1 if r['key'].startswith('FJ1486_') else 2 if r['key'].startswith('FJ1477_') else 0
    fid=sys.argv[sys.argv.index('--muscle')+1] if '--muscle' in sys.argv else 'FJ1486';designs={name:json.loads((DEST/name/'design.json').read_text(encoding='utf8')) for name in ['FJ1486','FJ1477']}
    def camera(name):
        if '--whole' in sys.argv:
            cam=mujoco.MjvCamera();cam.lookat[:]=[-.05,0,.02];cam.distance=1.1;cam.azimuth=-90 if name=='FJ1486' else 90;cam.elevation=-4;return cam
        p=designs[name]['selected'];cam=mujoco.MjvCamera();cam.lookat[:]=np.array(p['ball_centre_display_mm'])/1000;cam.distance=.115;u=np.array(p['neutral_terminal_axis']);cam.azimuth=np.degrees(np.arctan2(u[1],u[0]))+188;cam.elevation=-np.degrees(np.arcsin(u[2]))+5;return cam
    if '--animate' in sys.argv:
        from PIL import Image,ImageDraw
        from math import cos,sin,radians
        renderer=mujoco.Renderer(m,520,650);frames=[]
        sweep=list(range(26))+list(range(24,-1,-1))+list(range(-1,-26,-1))+list(range(-24,1))
        for direction in ['v','w']:
            for angle in sweep:
                d.mocap_pos[:]=0;d.mocap_quat[:]=[1,0,0,0]
                for name in ['FJ1486','FJ1477']:
                    p=designs[name]['selected'];C=np.array(p['ball_centre_display_mm'])/1000;u=np.array(p['neutral_terminal_axis']);v=np.array(p['retainer_v_axis']);axis=v if direction=='v' else np.cross(u,v);a=radians(angle)
                    K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]);R=np.eye(3)+sin(a)*K+(1-cos(a))*(K@K)
                    bid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,name+'_swivel_ball_crimp_blank_motion');mid=int(m.body_mocapid[bid]);d.mocap_pos[mid]=C-R@C;d.mocap_quat[mid]=np.r_[cos(a/2),axis*sin(a/2)]
                mujoco.mj_forward(m,d);panel=Image.new('RGB',(1300,520),'white')
                for i,name in enumerate(['FJ1486','FJ1477']):
                    option=mujoco.MjvOption();option.geomgroup[:]=0;option.geomgroup[i+1]=1;renderer.update_scene(d,camera=camera(name),scene_option=option);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
                    im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,650,60),fill='white');draw.text((10,10),'M5 '+name+' | prescribed local sweep '+str(angle)+'deg',fill='black');draw.text((10,32),'TPU fixed; no drive/force/deformation simulation',fill='black');panel.paste(im,(650*i,0))
                frames.append(panel.convert('P',palette=Image.Palette.ADAPTIVE,colors=160))
        output=DEST/'terminal_articulation.gif';frames[0].save(output,save_all=True,append_images=frames[1:],duration=85,loop=0,optimize=False);renderer.close()
        for filename,digest in inputs.items():assert hashlib.sha256(Path(filename).read_bytes()).hexdigest()==digest,'Rendering input changed'
        (DEST/'terminal_animation_inputs.json').write_text(json.dumps({'scope':'captured before model compilation, unchanged through prescribed animation','input_sha256':inputs,'GIF_sha256':hashlib.sha256(output.read_bytes()).hexdigest()},indent=2),encoding='utf8');print(output);return
    if '--render' in sys.argv:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(m,850,1000);panel=Image.new('RGB',(2000,850),'white')
        for i,name in enumerate(['FJ1486','FJ1477']):
            option=mujoco.MjvOption()
            if '--whole' not in sys.argv:option.geomgroup[:]=0;option.geomgroup[i+1]=1
            renderer.update_scene(d,camera=camera(name),scene_option=option);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1000,65),fill='white');version='M5B' if '--distal-load' in sys.argv else 'M5A';title=version+' | 131 actual neutral meshes' if '--whole' in sys.argv else version+' | '+name+' | isolated real terminal meshes';subtitle='Whole neutral packaging; no verified muscle motion or actuation' if '--whole' in sys.argv else 'Other anatomy hidden for detail; no TPU deformation or drive shown';draw.text((12,12),title,fill='black');draw.text((12,35),subtitle,fill='black');panel.paste(im,(1000*i,0))
        output=DEST/('whole_neutral.png' if '--whole' in sys.argv else 'terminal_cutaway.png' if '--cutaway' in sys.argv else 'terminal_detail.png');panel.save(output);renderer.close();print(output);return
    from mujoco import viewer
    import time
    cam=camera(fid)
    with viewer.launch_passive(m,d) as v:
        v.cam.lookat[:]=cam.lookat;v.cam.distance=cam.distance;v.cam.azimuth=cam.azimuth;v.cam.elevation=cam.elevation
        while v.is_running():v.sync();time.sleep(.03)
if __name__=='__main__':run()
