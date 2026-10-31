"""Chamfer interfering SC fork corners, with protected bearing carriers."""
import json,numpy as np
from verify_joint_candidate import load_candidate,OUT,world,md,mujoco
from build_candidate import inflated,to_local,sphere

def run():
    m,d,s,b=load_candidate();clav=s['clavicle_visual'];before=clav.volume()
    chest=inflated(s['thorax_geom'],.3)
    # Axis2 bearing carriers: local Z-axis, centers at +/-12 mm.
    protected=md.Manifold.batch_boolean([md.Manifold.cylinder(12,8,circular_segments=48,center=True).translate([0,0,z]) for z in [-12,12]],md.OpType.Add)
    collision_in_protected=0.
    knots=json.loads((OUT/'shoulder_path_fit.json').read_text())['knots'];qs=[np.zeros(3)]+[np.array(k['q_deg']) for k in knots if k['accepted_as_path_knot']]
    for q in qs:
        d.qpos[:]=0;d.qpos[:3]=np.radians(q);mujoco.mj_forward(m,d)
        cutter=to_local(d,world(d,chest,b['thorax_geom']),b['clavicle_visual'])^sphere(np.zeros(3),30)
        collision_in_protected=max(collision_in_protected,(clav^cutter^protected).volume())
        clav=clav-(cutter-protected)
    s['clavicle_visual']=clav
    rows=[];arrays={}
    for n,solid in s.items():
        simplified=solid.simplify(.02)
        rows.append({'bone':n,'triangles_before':solid.num_tri(),'triangles_after':simplified.num_tri(),'volume_before_simplification_mm3':solid.volume(),'volume_after_simplification_mm3':simplified.volume(),'positive_shell_volumes_mm3':[p.volume() for p in simplified.decompose() if p.volume()>0]})
        mesh=simplified.to_mesh();arrays[n]=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
    np.savez_compressed(OUT/'sc_repaired_local_meshes.npz',**arrays)
    report={'scope':'SC fork corner relief only; retained radius8 bearing zones around +/-12mm axis2 ears; surface simplification tolerance0.02mm must be followed by new full audit','clavicle_removed_volume_mm3':before-clav.volume(),'max_requested_cut_in_protected_carrier_mm3':collision_in_protected,'parts':rows,'frozen':False}
    (OUT/'sc_corner_repair.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':run()
