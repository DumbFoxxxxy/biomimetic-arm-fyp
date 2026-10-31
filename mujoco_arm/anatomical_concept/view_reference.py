"""Raw anatomical overlay, no kinematic model or soft simulation."""
from pathlib import Path
import argparse,time,json
import numpy as np,mujoco
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept'

def run():
    p=argparse.ArgumentParser()
    for flag in ['render','check','lattice','cage','attachments','supports','joints','SC','AC','GH','closed','compact','detail']:p.add_argument('--'+flag,action='store_true')
    p.add_argument('--preserve-head',action='store_true')
    p.add_argument('--elbow',action='store_true')
    p.add_argument('--muscle');args=p.parse_args();joint_view=args.joints or args.SC or args.AC or args.GH
    if args.elbow:args.GH=True;joint_view=True
    target=OUT/'joint_model/SC_connection/AC_connection' if args.AC else OUT/'joint_model/SC_connection' if args.SC else OUT/'joint_model' if args.joints else OUT/'attachment_regions' if args.attachments else OUT/'surface_cage_screen' if args.cage else OUT/'lattice_screen' if args.lattice else OUT
    if args.AC and args.compact:target=target.parent/'AC_compact'
    if args.GH:target=OUT/'joint_model/SC_connection/AC_compact'/('GH_closed' if args.closed else 'GH_preserve_head' if args.preserve_head else 'GH_connection')
    if args.elbow:target=target/'elbow_connection'
    m=mujoco.MjModel.from_xml_path(str(target/('arm.xml' if joint_view else 'reference.xml')));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    if args.SC or args.AC or args.GH:
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'FJ3290');m.geom_group[g]=3;m.geom_rgba[g,:]=[.5,.65,.7,1]
        if args.detail:
            for g in range(m.ngeom):
                name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g)
                keep=(name in ['FJ3384','FJ3368'] or name.startswith('GH_')) if args.GH else (name in ['FJ3362','FJ3384'] or name.startswith('AC_')) if args.AC else (name in ['FJ3362','FJ3290'] or name.startswith('SC_'))
                if args.elbow:keep=name in ['FJ3368','FJ3391'] or name.startswith('elbow_')
                if not keep:m.geom_rgba[g,3]=0
    if args.muscle:
        selected=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,args.muscle)
        if selected<0:raise ValueError('Unknown muscle '+args.muscle)
        for g in range(m.ngeom):
            if m.geom_group[g]==1 and g!=selected:m.geom_rgba[g,3]=0
        for s in range(m.nsite):
            if not mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_SITE,s).startswith(args.muscle+'_'):m.site_rgba[s,3]=0
    if args.attachments:
        m.geom_rgba[:,3]*=.2
        m.site_size[:]=.0025 # symbolic region samples, not physical anchor diameters
    if args.check:print('meshes',m.nmesh,'geoms',m.ngeom,'joints',m.njnt,'nv',m.nv);return
    def camera(cam):
        cam.lookat[:]=[-.015,-.005,.04];cam.distance=1.65 if args.supports else 1.35;cam.elevation=-4;cam.azimuth=-90
        if args.detail and (args.SC or args.AC or args.GH):
            fit=json.loads((OUT/'joint_candidates.json').read_text(encoding='utf8'));shift=np.array([-140.,-70.,1100.]);cam.lookat[:]=(np.array(fit['elbow' if args.elbow else 'GH' if args.GH else 'AC' if args.AC else 'SC']['centre_source_mm'])-shift)/1000;cam.distance=.3 if args.elbow else .24;cam.elevation=-15;cam.azimuth=-40
    if args.render:
        from PIL import Image,ImageDraw
        renderer=mujoco.Renderer(m,1000,1200);cam=mujoco.MjvCamera();camera(cam);views=[];options=mujoco.MjvOption()
        options.geomgroup[2]=int(args.supports)
        options.geomgroup[3]=int(args.SC or args.AC or args.GH)
        if joint_view:options.geomgroup[1]=0
        if args.attachments:options.sitegroup[3]=1
        for az,title in [(-90,'Anterior view'),(90,'Posterior view')]:
            cam.azimuth=(-40 if az==-90 else 140) if args.detail else az
            renderer.update_scene(d,camera=cam,scene_option=options);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,1200,64],fill='white')
            draw.text((20,12),'Independent anatomical reference | '+('Oblique detail '+str(az) if args.detail else title),fill='black')
            hardware='SC + AC + GH hardware' if args.GH else 'SC + AC hardware' if args.AC else 'SC hardware' if args.SC else 'ideal joints'
            if args.elbow:hardware='SC + AC + GH + elbow hardware'
            label='10DOF candidate; '+hardware+'; NO muscle drive validation' if joint_view else 'Atlas-constrained proximity regions, NOT verified attachments | '+str(args.muscle or 'all muscles') if args.attachments else 'Anatomical volume-clipped lattice screen | NOT fiber / wire / actuation validated' if args.cage or args.lattice else f'{sum(m.geom_group==0)} limb bones / {sum(m.geom_group==2) if args.supports else 0} support refs / {sum(m.geom_group==1)} muscles | anatomy only'
            draw.text((20,35),label,fill='black');im.save(target/('anatomy_'+str(az)+'.png'));views.append(im.resize((600,500)))
        panel=Image.new('RGB',(1200,500),'white');panel.paste(views[0],(0,0));panel.paste(views[1],(600,0));panel.save(target/'anatomical_overlay.png');renderer.close();print(target/'anatomical_overlay.png');return
    from mujoco import viewer as api
    play={'sequence':[],'index':0};review_paths={}
    if joint_view and (target/'bone_motion_audit.json').exists():
        from joint_motion_data import accepted_paths
        review_paths=accepted_paths(target)
    def key(k):
        if k in [ord('M'),ord('m')]:v.opt.geomgroup[1]^=1
        if k in [ord('B'),ord('b')]:v.opt.geomgroup[0]^=1
        if k in [ord('P'),ord('p')]:v.opt.sitegroup[3]^=1
        if k in [ord('T'),ord('t')]:v.opt.geomgroup[2]^=1
        if k in [ord('R'),ord('r')]:camera(v.cam);d.qpos[:]=m.qpos0;play['sequence']=[]
        if joint_view and k in [ord(str(i)) for i in range(1,6)]:
            name=['elbow_flexion','GH_flexion','GH_abduction','coordinated_raise','SC_elevation'][k-ord('1')]
            if name in review_paths:play['sequence']=review_paths[name];play['index']=0
    with api.launch_passive(m,d,key_callback=key) as v:
        camera(v.cam)
        v.opt.geomgroup[2]=int(args.supports)
        v.opt.geomgroup[3]=int(args.SC or args.AC or args.GH)
        if joint_view:v.opt.geomgroup[1]=0
        if args.attachments:v.opt.sitegroup[3]=1
        last=time.perf_counter()
        while v.is_running():
            if joint_view:
                seq=play['sequence']
                if seq and time.perf_counter()-last>.04:
                    i=min(play['index'],len(seq)-1);d.qpos[:]=seq[i];play['index']=(i+1)%len(seq);last=time.perf_counter()
                mujoco.mj_forward(m,d)
                if np.max(np.abs(d.qpos-m.qpos0))>1e-6:v.opt.geomgroup[1]=0
            v.sync();time.sleep(.03)

if __name__=='__main__':run()
