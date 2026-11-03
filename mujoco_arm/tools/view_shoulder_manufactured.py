"""Static viewer of the audited available parts, with recorded placements."""
from pathlib import Path
import argparse
import sys
import json
import xml.etree.ElementTree as ET
import struct
import hashlib
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from project_runtime import load_mujoco_model

DEST=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_integrated_revision'


def main():
    global DEST
    selection_file=ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection=json.loads(selection_file.read_text(encoding='utf8')) if selection_file.exists() else None
    selected_report=ROOT/selection['assembly_report'] if selection else DEST/'actual_assembly_check.json'
    ap=argparse.ArgumentParser()
    ap.add_argument('--report',type=Path,default=selected_report)
    ap.add_argument('--render',action='store_true');ap.add_argument('--check',action='store_true')
    ap.add_argument('--detail',action='store_true',help='Close-up of the FJ1467 outlet candidate')
    ap.add_argument('--part',help='Isolate one exact source key, for example cap:FJ1500/distal')
    ap.add_argument('--title',help='Explicit render label for a partial or comparison scene')
    ap.add_argument('--note',help='Explicit limitation text for the rendered scene')
    ap.add_argument('--azimuth',type=float,help='Optional inspection camera azimuth in degrees')
    ap.add_argument('--elevation',type=float,help='Optional inspection camera elevation in degrees')
    ap.add_argument('--distance',type=float,help='Optional inspection camera distance in metres')
    args=ap.parse_args();DEST=args.report.parent
    if selection and args.report.resolve()==selected_report.resolve() and hashlib.sha256(args.report.read_bytes()).hexdigest()!=selection['assembly_sha256']:
        raise ValueError('Current research-candidate selection is stale; update its verified report hashes')
    isolated=args.part or ('cap:FJ1467/origin' if args.detail else None)
    import mujoco
    import numpy as np
    data=json.loads(args.report.read_text(encoding='utf8'))
    root=ET.Element('mujoco',model='available_S1_parts_NOT_assembly_cleared')
    ET.SubElement(root,'compiler',angle='radian')
    visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='900')
    asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody')
    ET.SubElement(world,'light',pos='0 -1 2',dir='0 1 -1')
    for i,(key,row) in enumerate(data['sources'].items()):
        if isolated and key!=isolated:continue
        name='part_%03d'%i
        source=ROOT/row['file'];display=source
        actual_digest=hashlib.sha256(source.read_bytes()).hexdigest()
        if row.get('sha256') and actual_digest!=row['sha256']:
            raise ValueError('Changed actual scene input: '+key)
        if source.suffix.lower()=='.stl':
            with source.open('rb') as stream:
                stream.seek(80);count=struct.unpack('<I',stream.read(4))[0]
            if count>200000:
                # MuJoCo's STL decoder caps face count. Reindex the exact same
                # triangles as OBJ; never decimate a manufacturing part to make
                # the visual claim easier to render.
                folder=DEST/'view_assets';folder.mkdir(exist_ok=True)
                display=folder/(name+'.obj');record=folder/(name+'.json')
                digest=actual_digest
                if not display.exists() or not record.exists() or json.loads(record.read_text())['source_sha256']!=digest:
                    blob=source.read_bytes();dtype=np.dtype([('n','<f4',3),('v','<f4',(3,3)),('a','<u2')])
                    triangles=np.frombuffer(blob,dtype=dtype,offset=84)['v']
                    vertices,index=np.unique(triangles.reshape(-1,3),axis=0,return_inverse=True)
                    faces=index.reshape(-1,3)+1
                    with display.open('w',encoding='utf8') as stream:
                        for v in vertices:stream.write('v %.9g %.9g %.9g\n'%tuple(v))
                        for f in faces:stream.write('f %d %d %d\n'%tuple(f))
                    record.write_text(json.dumps(dict(source=str(source.relative_to(ROOT)),source_sha256=digest,triangles=count,vertices=len(vertices),decimated=False)),encoding='utf8')
        ET.SubElement(asset,'mesh',name=name,file=str(display),scale='.001 .001 .001')
        position=np.array(row['translation_source_mm'] or [0,0,0],float)*.001
        colour='.9 .87 .8 1' if key.startswith('bone:') else '.8 .12 .1 .85' if key.startswith('tpu:') else '.1 .7 .9 1' if key.startswith('tendon:') else '.6 .65 .7 1'
        ET.SubElement(world,'geom',name=name,type='mesh',mesh=name,rgba=colour,
                      pos=' '.join(map(str,position)),contype='0',conaffinity='0')
    detail_stem='part_detail' if args.part else 'outlet_detail'
    scene_file=DEST/(detail_stem+'.xml' if isolated else 'available_parts.xml')
    ET.ElementTree(root).write(scene_file,encoding='utf-8',xml_declaration=True)
    m=load_mujoco_model(scene_file);d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    if isolated and m.ngeom!=1:raise ValueError('Requested exact part not found: '+isolated)
    print('Actual available scene geoms',m.ngeom,'DOF',m.nq)
    def camera(cam):
        cam.lookat[:]=[-.14,-.04,1.285];cam.distance=.50;cam.azimuth=130;cam.elevation=6
        if isolated:cam.lookat[:]=d.geom_xpos[0];cam.distance=max(.035,float(m.geom_rbound[0])*6);cam.azimuth=90;cam.elevation=15
        if args.azimuth is not None:cam.azimuth=args.azimuth
        if args.elevation is not None:cam.elevation=args.elevation
        if args.distance is not None:cam.distance=args.distance
    if args.check:return
    if args.render:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(m,height=900,width=1200);cam=mujoco.MjvCamera();camera(cam)
        renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im)
        draw.rectangle((0,0,1200,65),fill='white')
        title=f'ISOLATED {isolated} | SURROUNDING ANATOMY HIDDEN' if isolated else f'ACTUAL AVAILABLE S1 PARTS | {m.ngeom} nominal instances | NOT ASSEMBLY-CLEARED'
        draw.text((12,10),args.title or title,fill='black')
        draw.text((12,35),args.note or 'Static manufacturing meshes; no validated motion, liner, fasteners or muscle deformation',fill='black')
        im.save(DEST/(detail_stem+'.png' if isolated else 'available_parts.png'));renderer.close();return
    import time
    from mujoco import viewer
    with viewer.launch_passive(m,d) as v:
        camera(v.cam)
        while v.is_running():v.sync();time.sleep(.03)


if __name__=='__main__':main()
