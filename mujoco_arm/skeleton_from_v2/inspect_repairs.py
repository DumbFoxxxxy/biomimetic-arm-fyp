"""Isolated local-shape experiments. Never edits V2 or declares a freeze."""
import json, sys
import numpy as np
from audit_baseline import ROOT, local_solid, mujoco
import manifold3d as md

NAMES=['thorax_geom','clavicle_visual','scapula_mesh','humerus_shaft','ulna_shaft']
def load():
    m=mujoco.MjModel.from_xml_path(str(ROOT/'mujoco_arm/cable_design_v2/arm_cable_design_v2.xml'))
    d=mujoco.MjData(m);mujoco.mj_forward(m,d)
    solids={};bodies={}
    for name in NAMES:
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,name)
        solids[name]=local_solid(m,d,g)[0];bodies[name]=int(m.geom_bodyid[g])
    return m,d,solids,bodies

def sphere(p,r):return md.Manifold.sphere(r,circular_segments=48).translate(p)

def shifted_plate(source,offset):
    # Normal of the brief's plane, pointing away from chest.
    ac=np.array([156.86,11.97,233.24]); gh=np.array([159.06,-7.30,187.54])
    sup=np.array([72.76,-82.21,268.45]);inf=np.array([44.55,-106.69,112.87])
    w=np.cross(inf-sup,gh-sup);w/=np.linalg.norm(w)
    if w[1]>0:w=-w
    centers=[np.zeros(3),gh-ac]
    def warp(p):
        p=np.array(p);dist=min(np.linalg.norm(p-c) for c in centers)
        t=np.clip((dist-18)/32,0,1);t=t*t*(3-2*t)
        return p+offset*w*t
    return source.warp(warp)

def world(d,solid,body):
    return solid.transform(np.column_stack([d.xmat[body].reshape(3,3),d.xpos[body]*1000]))

def run():
    m,d,original,bodies=load();rows=[]
    for offset in [0,4,8,12,16]:
        s=shifted_plate(original['scapula_mesh'],offset)
        overlaps={}
        for other in ['thorax_geom','clavicle_visual','humerus_shaft']:
            a=world(d,s,bodies['scapula_mesh']);b=world(d,original[other],bodies[other]);hit=a^b
            overlaps[other]={'volume_mm3':hit.volume(),'bounds':list(hit.bounding_box()) if not hit.is_empty() else None}
        row={'offset_mm':offset,'status':str(s.status()),'volume_mm3':s.volume(),'signed_shell_volumes':[x.volume() for x in s.decompose()],'zero_overlaps':overlaps}
        rows.append(row);print(json.dumps(row),flush=True)
    out=ROOT/'mujoco_arm/exports/skeleton_from_v2';out.mkdir(exist_ok=True,parents=True)
    (out/'local_repair_experiments.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':run()
