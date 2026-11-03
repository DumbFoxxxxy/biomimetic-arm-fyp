"""MuJoCo mocap display of independently checked bare GH geometry paths.

No muscles, collar, tendons, dynamic actuation or actual driven range.
Uses the exact manufacturing triangles via the existing static mesh loader.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from project_runtime import load_mujoco_model


def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def run(report_path,axis_label,sign,render=False,check=False):
    import numpy as np
    import mujoco
    report_path=report_path.resolve();report=json.loads(report_path.read_text(encoding='utf8'));out=report_path.parent
    source=ROOT/report['assembly_source']
    if hashed(source)!=report['assembly_sha256']:raise ValueError('Changed motion source assembly')
    data=json.loads(source.read_text(encoding='utf8'))
    path=next(r for r in report['results'] if r['axis']==axis_label and r['sign']==sign)
    maximum=path['last_contiguous_sample_clear_deg']
    if maximum is None or not path['return_sampled_path_clear']:raise ValueError('No checked outward/return path')
    keys=report['fixed_keys']+report['moving_keys']
    selected={k:data['sources'][k] for k in keys}
    if set(keys)&set(report['not_installed_keys']):raise ValueError('Stage membership mismatch')
    for k,r in selected.items():
        if hashed(ROOT/r['file'])!=r['sha256']:raise ValueError('Changed scene mesh '+k)
    scene=out/'bare_stage_sources.json';scene.write_text(json.dumps(dict(sources=selected),indent=2)+'\n',encoding='utf8')
    subprocess.run([sys.executable,str(ROOT/'mujoco_arm/tools/view_shoulder_manufactured.py'),
                    '--report',str(scene),'--check'],check=True)
    xml=ET.parse(out/'available_parts.xml');root=xml.getroot();root.set('model','bare_GH_mocap_NOT_actuation')
    ET.SubElement(root,'option',gravity='0 0 0')
    world=root.find('worldbody');index=keys.index('bone:humerus');name=f'part_{index:03d}'
    geom=world.find(f"geom[@name='{name}']")
    if geom is None:raise ValueError('Missing exact humerus geometry')
    centre=np.array(report['joint_centre_source_mm'])*.001;old=np.array([float(x) for x in geom.get('pos').split()])
    world.remove(geom);geom.set('pos',' '.join(map(str,old-centre)))
    body=ET.SubElement(world,'body',name='bare_humerus',mocap='true',pos=' '.join(map(str,centre)));body.append(geom)
    xml_path=out/'bare_GH_motion.xml';xml.write(xml_path,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(xml_path);state=mujoco.MjData(model);mujoco.mj_forward(model,state)
    if model.nq!=0 or model.nmocap!=1 or model.ngeom!=len(keys):raise ValueError('Unexpected display model scope')
    axis=np.array(report['axes_source'][axis_label]);axis/=np.linalg.norm(axis)
    gid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,name)
    initial=state.geom_xpos.copy();fixed_ids=[i for i in range(model.ngeom) if i!=gid]
    def pose(degrees):
        half=np.radians(degrees)/2.;q=np.r_[np.cos(half),axis*np.sin(half)]
        state.mocap_quat[0]=q;mujoco.mj_forward(model,state)
        a=np.radians(degrees);k=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0.]])
        R=np.eye(3)+np.sin(a)*k+(1-np.cos(a))*(k@k)
        expected=centre+R@(initial[gid]-centre)
        if not np.allclose(state.geom_xpos[gid],expected,atol=1e-9):raise ValueError('Display bone rotation mismatch')
        if not np.allclose(state.geom_xpos[fixed_ids],initial[fixed_ids],atol=1e-12):raise ValueError('Fixed scene moved')
    for value in [0.,maximum/2.,maximum]:pose(value)
    metadata=dict(scope=__doc__,motion_report_sha256=hashed(report_path),assembly_sha256=hashed(source),
                  script_sha256=hashed(Path(__file__).resolve()),axis=axis_label,sign=sign,
                  displayed_last_clear_deg=maximum,scene_geoms=model.ngeom,source_object_count=report['source_object_count'],
                  nq=model.nq,nmocap=model.nmocap,fixed_keys=report['fixed_keys'],not_installed_keys=report['not_installed_keys'],
                  geometry_rotation_checks_pass=True,rendered_frames=[0.,maximum/2.,maximum] if render else [],
                  stage_pass=False,actual_drive_range_deg=None)
    (out/'bare_GH_view_checks.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf8')
    print('Bare scene',model.ngeom,'geoms; mocap',model.nmocap,'DOF',model.nq,'angle',maximum,flush=True)
    if check:return
    def camera(cam):cam.lookat[:]=[-.14,-.04,1.285];cam.distance=.55;cam.azimuth=130;cam.elevation=6
    if render:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(model,height=700,width=900);cam=mujoco.MjvCamera();camera(cam)
        frames=[]
        for value in [0.,maximum/2.,maximum]:
            pose(value);renderer.update_scene(state,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
            im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,900,65),fill='white')
            draw.text((10,8),f'BARE GH {axis_label} {value:.1f} deg | prescribed geometry',fill='black')
            draw.text((10,32),'Muscles, collar and tendons absent; no actuation or strength acceptance',fill='black');frames.append(im)
        canvas=Image.new('RGB',(2700,700),'white')
        for i,im in enumerate(frames):canvas.paste(im,(i*900,0))
        canvas.save(out/'bare_GH_motion.png');renderer.close();return
    from mujoco import viewer
    playing=True;start=time.perf_counter();manual=None
    def keypress(key):
        nonlocal playing,manual,start
        if key==32:playing=not playing;manual=None;start=time.perf_counter()
        elif key==48:playing=False;manual=0.
        elif key==49:playing=False;manual=maximum
    with viewer.launch_passive(model,state,key_callback=keypress) as v:
        camera(v.cam)
        while v.is_running():
            angle=manual if manual is not None else maximum*(.5-.5*np.cos((time.perf_counter()-start)*np.pi/3)) if playing else 0.
            pose(angle);v.sync();time.sleep(.03)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--report',type=Path,required=True);parser.add_argument('--axis',choices=['flexion','abduction','longitudinal'],default='flexion')
    parser.add_argument('--sign',type=int,choices=[-1,1],default=1);parser.add_argument('--render',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();run(args.report,args.axis,args.sign,args.render,args.check)
