"""SC support swept slot with protected bearing zones and topology gate."""
import json,numpy as np
from verify_joint_candidate import load_candidate,OUT,world,md,mujoco
from build_candidate import inflated,to_local,sphere
def run():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--fit',default='shoulder_path_fit.json');args=parser.parse_args()
    m,d,s,b=load_candidate();clav=inflated(s['clavicle_visual'],.3);chest=s['thorax_geom'];before=chest.volume();sc=np.array([40,106,235])
    protected=md.Manifold.batch_boolean([md.Manifold.cylinder(8,8,circular_segments=48,center=True).rotate([90,0,0]).translate(sc+[0,y,0]) for y in [-12,12]],md.OpType.Add)
    knots=json.loads((OUT/args.fit).read_text())['knots'];protected_conflict=0
    for k in knots:
        if not k['accepted_as_path_knot']:continue
        d.qpos[:]=0;d.qpos[:3]=np.radians(k['q_deg']);mujoco.mj_forward(m,d)
        cutter=to_local(d,world(d,clav,b['clavicle_visual']),b['thorax_geom'])^sphere(sc,32)
        protected_conflict=max(protected_conflict,(chest^cutter^protected).volume())
        chest=chest-(cutter-protected)
    s['thorax_geom']=chest;arrays={};parts=[]
    for n,solid in s.items():
        simplified=solid.simplify(.02);components=[x.volume() for x in simplified.decompose() if x.volume()>1]
        parts.append({'bone':n,'volume_mm3':simplified.volume(),'positive_components_over_1mm3':components,'triangles':simplified.num_tri(),'status':str(simplified.status())})
        if len(components)!=1:raise RuntimeError('Disconnected structural components: '+n)
        mesh=simplified.to_mesh();arrays[n]=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
    np.savez_compressed(OUT/'sc_mount_repaired_local_meshes.npz',**arrays)
    report={'frozen':False,'fit_file':args.fit,'SC_center_moved_mm':0,'chest_local_removed_volume_mm3':before-chest.volume(),'max_protected_bearing_conflict_mm3':protected_conflict,'surface_simplification_tolerance_mm':.02,'parts':parts}
    (OUT/'sc_mount_repair.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':run()
