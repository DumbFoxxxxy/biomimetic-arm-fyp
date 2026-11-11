"""Render zero-DOF rest reference, including the unchanged original elbow lattices."""
import sys
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from render_channel_first_conduits import png


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_elbow_supply'
    scene_file=folder/'neutral_scene.json';scene=read(scene_file)
    model_file=ROOT/scene['model']['file'];assert provenance(model_file)==scene['model']
    model=load_mujoco_model(model_file);assert model.nv==0
    data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    renderer=mujoco.Renderer(model,height=900,width=1100);outputs=[]
    for label,target,distance in [('neutral_whole',[-.08,-.13,1.12],.90),('neutral_routes',[-.12,-.13,1.27],.50)]:
        camera=mujoco.MjvCamera();camera.lookat[:]=target;camera.distance=distance;camera.azimuth=55.;camera.elevation=8.
        renderer.update_scene(data,camera=camera);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False
        output=folder/(label+'.png');png(output,renderer.render());outputs.append(provenance(output))
    renderer.close();save(folder/'render_record.json',dict(scene=provenance(scene_file),images=outputs,
        generator=provenance(ROOT/'tools/render_channel_first_external_supply.py'),
        scope='Actual MuJoCo neutral geometry rendering: green unloaded hollow supply references, red original elbow lattice meshes unchanged, blue local sternum frame. ZeroDOF; no actuator, deformation or dynamic supply approval.'))
    print('Neutral whole-arm and supply reference views rendered',flush=True)


if __name__=='__main__':main()
