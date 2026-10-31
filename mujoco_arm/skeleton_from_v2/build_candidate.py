"""Local relief candidate, retained V2 centers/chain. Not a frozen design."""
import itertools,json
import numpy as np
from inspect_repairs import load,world,sphere,shifted_plate,ROOT,md,mujoco,NAMES
from pathlib import Path
import sys
sys.path.insert(0,str(ROOT/'mujoco_arm'))
from stl_io import write_binary_stl
OUT=ROOT/'mujoco_arm/exports/skeleton_from_v2'

def inflated(s,margin=.3):
    # Conservative box dilation. Exact convex hull of translated copies is
    # inappropriate for concave bones; use Minkowski sum with a small cube.
    return s.minkowski_sum(md.Manifold.cube([2*margin]*3,center=True))

def to_local(d,solid,b):
    r=d.xmat[b].reshape(3,3);return solid.transform(np.column_stack([r.T,-r.T@d.xpos[b]*1000]))

def build(offset=8):
    m,d,s,b=load();before={n:x for n,x in s.items()}
    # Refine planar source faces before smooth deformation so a coarse long
    # triangle cannot turn stand-off into spurious plate thickening.
    s['scapula_mesh']=shifted_plate(s['scapula_mesh'].refine_to_length(4),offset)
    ac=d.xpos[b['scapula_mesh']]*1000;gh=d.xpos[b['humerus_shaft']]*1000;e=d.xpos[b['ulna_shaft']]*1000
    # Bring the AC ball up to the brief's R10 nominal, retaining its center.
    s['clavicle_visual']=s['clavicle_visual']+to_local(d,sphere(ac,10),b['clavicle_visual'])
    scap=world(d,s['scapula_mesh'],b['scapula_mesh'])
    # True cavities/reliefs, not collision exclusions. Only local joint zones
    # are cut. Their topology and retained material are recorded for review.
    clav=world(d,s['clavicle_visual'],b['clavicle_visual'])
    hum=world(d,s['humerus_shaft'],b['humerus_shaft'])
    scap=scap-inflated(world(d,s['thorax_geom'],b['thorax_geom']))
    scap=scap-(inflated(clav)^sphere(ac,22))
    scap=scap-(inflated(hum)^sphere(gh,34))
    # Remove only isolated machining chips under 1 mm3, never a structural
    # joint component. A larger separated region makes this candidate fail.
    shells=scap.decompose()
    chips=[x for x in shells if 0<x.volume()<1]
    if chips:scap=scap-md.Manifold.batch_boolean(chips,md.OpType.Add)
    s['scapula_mesh']=to_local(d,scap,b['scapula_mesh'])
    ulna=world(d,s['ulna_shaft'],b['ulna_shaft'])
    ulna=ulna-(inflated(hum)^sphere(e,28))
    s['ulna_shaft']=to_local(d,ulna,b['ulna_shaft'])
    return m,d,s,b,before

def run():
    m,d,s,b,before=build();OUT.mkdir(exist_ok=True,parents=True);rows=[];arrays={}
    for name,solid in s.items():
        mesh=solid.to_mesh();tri=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
        write_binary_stl(OUT/(name+'_candidate_local.stl'),tri*.001,'V2 local repair candidate, body local metres')
        arrays[name]=np.asarray(world(d,solid,b[name]).to_mesh().vert_properties)[:,:3][np.asarray(world(d,solid,b[name]).to_mesh().tri_verts)]
        rows.append({'bone':name,'status':str(solid.status()),'before_volume_mm3':before[name].volume(),'after_volume_mm3':solid.volume(),'signed_shell_volumes_mm3':[x.volume() for x in solid.decompose()]})
    np.savez_compressed(OUT/'candidate_world_meshes.npz',**arrays)
    report={'candidate_only':True,'frozen':False,'offset_mm':8,'joint_centers_moved_mm':0,'changes':rows,'zero_overlaps':[]}
    for a,c in itertools.combinations(NAMES,2):
        hit=world(d,s[a],b[a])^world(d,s[c],b[c]);v=hit.volume()
        if v>.001:report['zero_overlaps'].append({'bones':[a,c],'volume_mm3':v})
    (OUT/'candidate_geometry.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
if __name__=='__main__':run()
