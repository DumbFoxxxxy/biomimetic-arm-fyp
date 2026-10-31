"""Budgeted reverse-motion contour relief with protected joint carriers."""
import json,itertools,numpy as np
from verify_joint_candidate import load_candidate,OUT,world,NAMES,md,mujoco
from build_candidate import sphere,inflated,to_local
def run():
    m,d,s,b=load_candidate('sc_mount_repaired_local_meshes.npz');scap=s['scapula_mesh'];before=scap.volume()
    chest_source=s['thorax_geom'];sc=np.array([40,106,235]);bearing=md.Manifold.batch_boolean([md.Manifold.cylinder(8,8,circular_segments=48,center=True).rotate([90,0,0]).translate(sc+[0,y,0]) for y in [-12,12]],md.OpType.Add)
    clav_cutter=inflated(s['clavicle_visual'],.3)
    for angle in np.arange(-.5,-6.1,-.5):
        d.qpos[:]=0;d.qpos[1]=np.radians(angle);mujoco.mj_forward(m,d)
        cutter=to_local(d,world(d,clav_cutter,b['clavicle_visual']),b['thorax_geom'])^sphere(sc,32)
        trial=s['thorax_geom']-(cutter-bearing)
        if len([p for p in trial.decompose() if p.volume()>1])!=1:break
        if chest_source.volume()-trial.volume()>1000:break
        s['thorax_geom']=trial
    chest=inflated(s['thorax_geom'],.3)
    protect=sphere(np.zeros(3),12)+sphere(np.array([2.2,-19.27,-45.7]),14);rows=[]
    floor=.85*83684.36577568574
    for name,axis,end in [('SLIDE_BACK',1,-6),('SCAP_INT',2,-20)]:
        for angle in np.arange(-.5,end-.1,-.5):
            q=np.zeros(7);q[axis]=angle;d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d)
            cutter=to_local(d,world(d,chest,b['thorax_geom']),b['scapula_mesh'])
            trial=scap-(cutter-protect);positive=[p.volume() for p in trial.decompose() if p.volume()>1]
            reason=None
            if trial.volume()<floor:reason='15 percent total scapular material-removal budget'
            elif len(positive)!=1:reason='structural material would disconnect'
            pairs=[]
            for a,c in itertools.combinations(NAMES,2):
                sa=trial if a=='scapula_mesh' else s[a];sb=trial if c=='scapula_mesh' else s[c]
                hit=world(d,sa,b[a])^world(d,sb,b[c])
                if hit.volume()>.001:pairs.append({'bones':[a,c],'volume_mm3':hit.volume()})
            if pairs:reason=reason or 'remaining collision outside allowed contour relief'
            row={'action':name,'trial_joint_deg':float(angle),'scapular_volume_mm3':trial.volume(),'accepted':reason is None,'reason':reason,'remaining_collisions':pairs}
            rows.append(row)
            if reason:break
            scap=trial;s['scapula_mesh']=scap
    s['scapula_mesh']=scap;arrays={}
    for n,solid in s.items():
        if n=='scapula_mesh':solid=solid.simplify(.01)
        mesh=solid.to_mesh();arrays[n]=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
    np.savez_compressed(OUT/'reverse_refined_local_meshes.npz',**arrays)
    report={'frozen':False,'additional_SC_mount_reverse_slot_mm3':chest_source.volume()-s['thorax_geom'].volume(),'additional_scapular_removal_mm3':before-scap.volume(),'total_scapular_removed_fraction':1-scap.volume()/83684.36577568574,'protected_AC_radius_mm':12,'protected_GH_radius_mm':14,'trials':rows}
    (OUT/'reverse_clearance_repair.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':run()
