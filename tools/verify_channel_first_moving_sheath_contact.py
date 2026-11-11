"""Independent world-triangle and MuJoCo applyFT checks at a nonzero test pose."""
import sys
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save,load_component,mesh_geom
from build_channel_first_inventory import provenance
from channel_first_moving_sheath_contact import MovingSheathContact


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_elbow_supply_dynamic/axial15_bend15'
    check_file=folder/'model_check.json';check=read(check_file);response_file=folder/'short_response.json';response=read(response_file)
    model=load_mujoco_model(ROOT/check['model']['file']);data=mujoco.MjData(model)
    state=min(response['trace'],key=lambda r:abs(r['time_s']-.04002))
    mujoco.mj_integratePos(model,data.qpos,np.array(state['generalized_displacement']),1.)
    # Synthetic small velocities exercise friction signs; this is verification,
    #not another physical trajectory or a prescribed-motion drive result.
    data.qvel[:]=np.random.default_rng(8419).normal(0.,.001,model.nv);mujoco.mj_forward(model,data)
    contact=MovingSheathContact(model,check,friction=.2);data.qfrc_applied[:]=0.;contact.apply(data)
    assert contact.latest,'Test pose did not engage actual surface contact'
    expected=np.zeros(model.nv);distance_error=0.;closest_error=0.;rows=[]
    beams={b['head']:b for b in contact.beams};bones={b['key']:b for b in contact.bones}
    world_triangles={}
    for item in contact.latest:
        beam=beams[item['beam']];bone=bones[item['bone']];i=item['station'];f=float(beam['fractions'][i])
        sa=beam['ids'][beam['a'][i]];sb=beam['ids'][beam['b'][i]];a=data.site_xpos[sa];b=data.site_xpos[sb]
        centre=(1-f)*a+f*b;point=np.array(item['point_source_mm'])*.001;force=np.array(item['force_on_guide_N'])
        torque=np.cross(point-centre,force)
        for weight,site in [(1-f,sa),(f,sb)]:
            mujoco.mj_applyFT(model,data,force*weight,torque*weight,data.site_xpos[site],int(model.site_bodyid[site]),expected)
        mujoco.mj_applyFT(model,data,-force,np.zeros(3),point,bone['body'],expected)
        if item['bone'] not in world_triangles:
            R=data.xmat[bone['body']].reshape(3,3);translation=data.xpos[bone['body']]*1000-R@bone['origin_mm']
            solid=bone['surface'].solid.transform(np.c_[R,translation]);vertices,faces=mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())
            world_triangles[item['bone']]=vertices[faces]
        tri=world_triangles[item['bone']];distances,closest=mesh_geom._closest_all(centre*1000,tri);j=int(np.argmin(distances))
        brute=float(distances[j]);reported=1.5-item['penetration_mm']
        de=abs(brute-abs(reported));qe=float(np.linalg.norm(closest[j]-point*1000))
        distance_error=max(distance_error,de);closest_error=max(closest_error,qe)
        rows.append(dict(bone=item['bone'],station=i,world_triangle_distance_error_mm=de,
            closest_point_error_mm=qe,tangential_power_W=item['tangential_power_W']))
    force_error=float(np.max(abs(expected-contact.held),initial=0))
    assert force_error<1e-8 and distance_error<1e-6 and closest_error<1e-6
    assert all(r['tangential_power_W']<=1e-12 for r in rows)
    save(folder/'moving_contact_verification.json',dict(model_check=provenance(check_file),pose_source=provenance(response_file),
        test_pose_time_s=state['time_s'],synthetic_velocity_seed=8419,contacts=rows,
        independent_applyFT_error_max_abs_mixed_units=force_error,
        brute_world_distance_error_max_mm=distance_error,brute_world_closest_point_error_max_mm=closest_error,
        generator=provenance(ROOT/'tools/verify_channel_first_moving_sheath_contact.py'),
        contact_module=provenance(ROOT/'tools/channel_first_moving_sheath_contact.py'),
        scope='Checks moving source/world transform against brute triangles and force/moment mapping against independent MuJoCo applyFT. Synthetic velocities only for friction dissipation. Does not validate material/contact stiffness, trajectory or continuum deformation.'))
    print('Moving contact verification',len(rows),'contacts;applyFT error',force_error,'world distance error mm',distance_error,flush=True)


if __name__=='__main__':main()
