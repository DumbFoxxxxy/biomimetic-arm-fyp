"""Actual tendon and four eye parts from a checked diagnostic route scene.

Anatomy is hidden. Prescribed static eye angles, no stepping, dynamics,
TPU deformation, cable tension or whole-shoulder movement is simulated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from project_runtime import load_mujoco_model
BASE=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/free_strand_routes'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--head',choices=['FJ1500','FJ1504'],default='FJ1504')
    parser.add_argument('--case',choices=['neutral','prescribed_eye_witness'],default='prescribed_eye_witness')
    parser.add_argument('--check',action='store_true');parser.add_argument('--render',action='store_true');args=parser.parse_args()
    report=json.loads((BASE/'route_audit.json').read_text(encoding='utf8'))
    checked=json.loads((BASE/'verification_checks.json').read_text(encoding='utf8'))
    if checked['audit_sha256']!=sha(BASE/'route_audit.json'):raise ValueError('Stale route verification')
    row=next(r for h in report['heads'] if h['head']==args.head for r in h['diagnostic_exports'] if r['label']==args.case)
    source=ROOT/row['scene'];scene=json.loads(source.read_text(encoding='utf8'))
    verification=next(r for r in checked['diagnostic_export_checks'] if r['head']==args.head and r['label']==args.case)
    if verification['scene_sha256']!=sha(source):raise ValueError('Changed diagnostic scene')
    import mujoco
    import numpy as np
    xml=ET.Element('mujoco',model='prescribed_static_tendon_witness_anatomy_hidden');ET.SubElement(xml,'compiler',angle='radian')
    visual=ET.SubElement(xml,'visual');ET.SubElement(visual,'global',offwidth='1000',offheight='800')
    asset=ET.SubElement(xml,'asset');world=ET.SubElement(xml,'worldbody');ET.SubElement(world,'light',pos='-.1 -.4 1.5',dir='0 1 -1')
    keys=list(row['exported_sources']);wire_name=None
    for i,key in enumerate(keys):
        entry=scene['sources'][key];file=ROOT/entry['file']
        if sha(file)!=entry['sha256']:raise ValueError('Changed displayed source')
        name='part_'+str(i);ET.SubElement(asset,'mesh',name=name,file=str(file),scale='.001 .001 .001')
        wire=key.endswith('/tendon_wire')
        if wire:wire_name=name
        colour='.95 .6 .12 1' if wire else '.12 .52 .74 1' if 'groove_ring' in key else '.65 .7 .74 1'
        ET.SubElement(world,'geom',name=name,type='mesh',mesh=name,rgba=colour,
            pos=' '.join(str(float(v)*.001) for v in entry['translation_source_mm']),contype='0',conaffinity='0')
    xmlfile=source.parent/'tendon_detail.xml';ET.ElementTree(xml).write(xmlfile,encoding='utf8',xml_declaration=True)
    model=load_mujoco_model(xmlfile);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    if model.ngeom!=5 or model.nq!=0:raise ValueError('Unexpected diagnostic scene')
    wire_id=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,wire_name)
    def camera(cam):
        cam.lookat[:]=data.geom_xpos[wire_id];cam.distance=max(.13,float(model.geom_rbound[wire_id])*5.5)
        cam.azimuth=110;cam.elevation=15
    print('Static route parts',model.ngeom,'DOF',model.nq,'eye angles',scene['phiA_deg'],scene['phiB_deg'])
    if args.check:return
    if args.render:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(model,height=800,width=1000);cam=mujoco.MjvCamera();camera(cam)
        renderer.update_scene(data,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1000,80),fill='white')
        draw.text((12,10),f'{args.head} | FULL TENDON + 4 EYE PARTS | ANATOMY HIDDEN',fill='black')
        draw.text((12,33),f'Prescribed eye angles {scene["phiA_deg"]:g}, {scene["phiB_deg"]:g} deg. Nominal length assumption.',fill='black')
        draw.text((12,56),'Static diagnostic geometry; no cable tension, TPU/cap equilibrium or shoulder drive result.',fill='black')
        im.save(source.parent/'tendon_detail.png');renderer.close();return
    from mujoco import viewer
    import time
    with viewer.launch_passive(model,data) as view:
        camera(view.cam)
        while view.is_running():view.sync();time.sleep(.03)

if __name__=='__main__':main()
