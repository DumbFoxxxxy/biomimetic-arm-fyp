"""Render the actual local sternum frame without a backboard or pretend actuators."""
import sys
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from render_channel_first_conduits import png


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_reaction_mount'
    check=read(folder/'model_check.json');model_file=ROOT/check['model']['file']
    assert provenance(model_file)==check['model']
    model=load_mujoco_model(model_file);model.vis.global_.offwidth=1100;model.vis.global_.offheight=900
    for row in check['geom_owners']:
        i=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,row['geom']);key=row['source_key']
        if key.startswith(('thorax_mount:','thorax_port:')):model.geom_rgba[i]=[.1,.55,.95,1.]
        elif key.startswith(('thorax_reference:','bone:')):model.geom_rgba[i]=[.8,.82,.75,.65]
    data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    renderer=mujoco.Renderer(model,height=900,width=1100);images=[]
    for label,target,distance in [('frame_close',[-.025,-.22,1.25],.30),('whole_arm',[-.08,-.08,1.15],.9)]:
        cam=mujoco.MjvCamera();cam.lookat[:]=target;cam.distance=distance;cam.azimuth=55.;cam.elevation=8.
        renderer.update_scene(data,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        output=folder/(label+'.png');png(output,renderer.render());images.append(provenance(output))
    renderer.close();save(folder/'render_record.json',dict(model=provenance(model_file),images=images,
        generator=provenance(ROOT/'tools/render_channel_first_thorax_mount.py'),
        scope='Actual235 rigid meshes, local frame blue, no backboard or motor. Green stubs neutral references; no complete supply or full TPU muscles shown.'))
    print('Sternum frame and whole arm rendered',flush=True)


if __name__=='__main__':main()
