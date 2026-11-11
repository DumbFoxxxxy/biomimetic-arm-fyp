"""Render actual inlet hardware and nominal sheath references from MuJoCo."""
import sys
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from render_channel_first_conduits import png


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface'
    check=read(folder/'model_check.json');assembly=read(folder/'assembly.json')
    model_file=ROOT/check['model']['file'];assert provenance(model_file)==check['model']
    model=load_mujoco_model(model_file);model.vis.global_.offwidth=960;model.vis.global_.offheight=720
    for row in check['geom_owners']:
        i=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,row['geom']);k=row['source_key']
        if k.endswith('proximal_rigid_cap') or k.endswith('_retainer'):model.geom_rgba[i]=[1.,.35,.05,1.]
        elif k.startswith('supply:'):model.geom_rgba[i]=[.72,.78,.83,1.]
        elif k.startswith(('bone:','thorax_reference:')):model.geom_rgba[i]=[.72,.76,.7,.4]
    data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    renderer=mujoco.Renderer(model,height=720,width=960);outputs=[]
    for row in assembly['interfaces']:
        target=np.array(row['origin_source_mm'])+np.array(row['reaction_axis_source'])*25.
        for azimuth in (55.,145.):
            cam=mujoco.MjvCamera();cam.lookat[:]=target*.001
            cam.distance=.12;cam.azimuth=azimuth;cam.elevation=12.
            renderer.update_scene(data,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
            image_file=folder/f"{row['id']}_inlet_{int(azimuth)}.png";png(image_file,renderer.render())
            outputs.append(dict(image=provenance(image_file),id=row['id'],azimuth_deg=azimuth,
                                lookat_m=cam.lookat.tolist(),distance_m=cam.distance))
    renderer.close()
    save(folder/'render_record.json',dict(model=provenance(model_file),assembly=provenance(folder/'assembly.json'),
        generator=provenance(ROOT/'tools/render_channel_first_elbow_supply.py'),images=outputs,
        scope='Neutral rendering of actual208rigid parts plus display-only short sheath references in green. Full chest routes and TPU muscle bodies absent. No deformation or force acceptance.'))
    print('Four actual inlet views rendered',flush=True)


if __name__=='__main__':main()
