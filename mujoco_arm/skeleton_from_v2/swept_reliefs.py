"""Evaluate physical swept reliefs, retaining material/connectivity diagnostics.

These are geometric experiments, not permission to hollow away arbitrary bones.
Final acceptance must inspect topology, shape loss, and bearing interfaces.
"""
import json,itertools
import numpy as np
from build_candidate import *

def apply_reliefs(m,d,s,b,scap_q,chest_sweep=True):
    original={k:v for k,v in s.items()};gh=np.array([2.2,-19.27,-45.7]);e=np.array([25.28,0,-288.89])
    def pose(q):d.qpos[:]=np.radians(q);mujoco.mj_forward(m,d)
    def relative(part,carrier):
        return to_local(d,world(d,s[part],b[part]),b[carrier])
    # GH/AC neck reliefs are limited to the immediate articulation region.
    cutters={n:inflated(s[n],.8 if n=='humerus_shaft' else .5) for n in ['humerus_shaft','clavicle_visual','ulna_shaft','thorax_geom']}
    def relcut(part,carrier):return to_local(d,world(d,cutters[part],b[part]),b[carrier])
    scap=s['scapula_mesh']
    for axis,end in [(3,120),(4,120),(5,70),(5,-90)]:
        for angle in np.linspace(0,end,int(abs(end)/2)+1):
            q=np.zeros(7);q[axis]=angle;pose(q)
            scap=scap-(relcut('humerus_shaft','scapula_mesh')^sphere(gh,39))
        print('GH relief',axis,end,scap.volume(),flush=True)
    # 2-degree swept cutters dilated by 0.8 mm cover the local interpolation
    # chord gap at radius38 (<0.67mm); all final paths are still audited.
    s['scapula_mesh']=scap
    cuts=[]
    for angle in np.linspace(-20,20,41):
        pose([0,0,angle,0,0,0,0]);s['scapula_mesh']=s['scapula_mesh']-(relcut('clavicle_visual','scapula_mesh')^sphere(np.zeros(3),23))
    print('AC relief',s['scapula_mesh'].volume(),flush=True)
    # Elbow fork/neck interference; retain original ulna and modify only the
    # humeral distal 55mm region. A disconnected fork is rejected downstream.
    cuts=[]
    for angle in np.linspace(0,145,146):
        pose([0,0,0,0,0,0,angle]);s['humerus_shaft']=s['humerus_shaft']-(relcut('ulna_shaft','humerus_shaft')^sphere(e,56))
    print('Elbow relief',s['humerus_shaft'].volume(),flush=True)
    # Shoulder paths: literal chest avoidance, not an ignored collision pair.
    cuts=[]
    targets=[np.r_[scap_q,np.zeros(4)],np.array([0,5,0,0,0,0,0]),np.array([0,-5,0,0,0,0,0]),np.array([0,0,20,0,0,0,0]),np.array([0,0,-20,0,0,0,0])]
    for q in targets if chest_sweep else []:
        for f in np.linspace(0,1,int(np.ceil(np.max(np.abs(q))))+1):
            pose(q*f);s['scapula_mesh']=s['scapula_mesh']-relcut('thorax_geom','scapula_mesh')
        print('Chest relief',q.tolist(),s['scapula_mesh'].volume(),flush=True)
    diagnostics=[]
    for name in ['scapula_mesh','humerus_shaft']:
        shell=s[name].decompose();positive=[p for p in shell if p.volume()>1]
        diagnostics.append({'bone':name,'volume_before_mm3':original[name].volume(),'volume_after_mm3':s[name].volume(),'removed_fraction':1-s[name].volume()/original[name].volume(),'positive_shells_over_1mm3':[p.volume() for p in positive],'status':str(s[name].status())})
        chips=[p for p in shell if 0<p.volume()<=1]
        if chips:s[name]=s[name]-md.Manifold.batch_boolean(chips,md.OpType.Add)
    pose(np.zeros(7));return diagnostics

if __name__=='__main__':
    m,d,s,b,_=build();rom=json.loads((OUT/'candidate_rom.json').read_text());q=rom['scapula_endpoint_fit']['q_deg']
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--joint-only',action='store_true');args=parser.parse_args()
    report=apply_reliefs(m,d,s,b,q,chest_sweep=not args.joint_only);print(json.dumps(report),flush=True)
    prefix='joint_relief' if args.joint_only else 'swept'
    (OUT/(prefix+'_experiment.json')).write_text(json.dumps(report,indent=2))
    arrays={}
    for name,solid in s.items():
        mesh=solid.to_mesh();tri=np.asarray(mesh.vert_properties)[:,:3][np.asarray(mesh.tri_verts)]
        arrays[name]=tri
        write_binary_stl(OUT/(name+'_'+prefix+'_local.stl'),tri*.001,'V2 swept-relief experiment, body local metres')
    np.savez_compressed(OUT/(prefix+'_local_meshes.npz'),**arrays)
