"""Actual MuJoCo packing render; no fictitious dynamic muscle constraints."""
from pathlib import Path
import json, hashlib, xml.etree.ElementTree as ET, sys
import numpy as np
import mujoco
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'mujoco_arm/exports/anatomical_concept';DEST=OUT/'elbow_guides'
CURRENT=OUT/'elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15/swivel_terminal_distal_load'

def run():
    grounded='--grounded' in sys.argv
    output=DEST/('ground_frame_swept' if '--swept-root' in sys.argv else 'ground_frame_low_root' if '--low-root' in sys.argv else 'ground_frame_candidate') if grounded else DEST
    selected=json.loads((DEST/'paired_placement/selection.json').read_text())['placements'];assembly_path=(output if grounded else CURRENT)/'assembly.json';assembly=json.loads(assembly_path.read_text());root=ET.Element('mujoco',model='M7 grounded guide packing candidate' if grounded else 'M6 ungrounded guide packing')
    ET.SubElement(root,'compiler',angle='radian');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1400',offheight='900');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 2',dir='0 1 -2',diffuse='.85 .85 .85');hashes={}
    def add(key,path,pos=(0,0,0),rgba='.78 .8 .82 1',body=None):
        path=path.resolve();hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest();ET.SubElement(asset,'mesh',name=key,file=str(path),scale='.001 .001 .001');ET.SubElement(body if body is not None else world,'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,pos)),rgba=rgba,contype='0',conaffinity='0')
    for r in assembly['instances']:
        body=None;key=r['key'];colour='.8 .81 .76 1' if key.startswith('FJ') and r['owner'] not in ['soft_unvalidated','distal_unvalidated'] else '.65 .68 .7 1'
        if r['kind']=='TPU':colour='.75 .1 .06 1'
        for fid,s in selected.items():
            if key!=fid+'_swivel_ball_crimp_blank':continue
            d=json.loads((CURRENT/fid/'design.json').read_text());C=np.array(d['selected']['ball_centre_display_mm']);u=np.array(d['selected']['neutral_terminal_axis']);L=s['source_distance_mm'];b=np.radians(s['offset_angle_deg']);p=-L*np.array([np.sin(b),np.cos(b)]);radius=8.005;t=radius**2/L**2*p-radius*np.sqrt(L**2-radius**2)/L**2*np.array([-p[1],p[0]]);incoming=(t-p)/np.linalg.norm(t-p);direction=np.array(s['frame_columns'])@np.r_[incoming,0];k=np.cross(u,direction);K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]]);R=np.eye(3)+K+K@K/(1+u@direction);q=np.zeros(4);mujoco.mju_mat2Quat(q,R.reshape(-1));body=ET.SubElement(world,'body',name=fid+'_terminal_alignment',pos=' '.join(map(str,(C-R@C)/1000)),quat=' '.join(map(str,q)))
        add(key,ROOT/r['file'],np.array(r['translation_display_mm'])/1000,colour,body)
    for path in sorted((DEST/'paired_placement').glob('*.obj')):
        if grounded and not path.stem.endswith('wire_corridor'):continue
        key=path.stem;colour='1 .25 .02 1' if key.endswith('wire_corridor') else '.1 .14 .19 1' if key.endswith('grooved_wheel') else '.25 .55 .85 1';add(key,path,rgba=colour)
    xml=output/'whole_arm.xml' if grounded else DEST/'paired_placement/whole_arm.xml';ET.ElementTree(root).write(xml,encoding='utf8',xml_declaration=True);m=mujoco.MjModel.from_xml_path(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    renderer=mujoco.Renderer(m,height=900,width=1400)
    for name,lookat,distance in [('whole_arm',[-.05,0,.025],1.05),('elbow_detail',[-.07,0,-.09],.43)]:
        cam=mujoco.MjvCamera();cam.lookat[:]=lookat;cam.distance=distance;cam.azimuth=125;cam.elevation=-7;renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1400,64),fill='white');draw.text((12,10),'M7 | real adapter-integrated guide frame candidate | MuJoCo' if grounded else 'M6 | actual arm + two sheave placement candidates | MuJoCo',fill='black');draw.text((12,34),'Ground-frame geometry candidate; cable ends FREE; TPU fixed; NO drive simulation' if grounded else 'Guides UNGROUNDED; cable ends FREE; TPU fixed; no contraction or drive simulated',fill='black');im.save(output/(name+'.png'))
    renderer.close()
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    (output/'render_inputs.json').write_text(json.dumps({'scope':'actual mesh inputs unchanged through static rendering','meshes':len(hashes),'input_sha256':hashes,'XML_sha256':hashlib.sha256(xml.read_bytes()).hexdigest(),'assembly_sha256':hashlib.sha256(assembly_path.read_bytes()).hexdigest()},indent=2),encoding='utf8');print('rendered real meshes',len(hashes))

if __name__=='__main__':run()
