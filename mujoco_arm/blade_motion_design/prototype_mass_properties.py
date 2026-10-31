"""Compute a consistent approximate printed-part inertia from each material union."""
import json,xml.etree.ElementTree as ET
import numpy as np
import mujoco
from shoulder_model import v8,EXPORT

def assign(root,density=1100.):
    from audit_gimbal_assembly import local_solid,md
    m=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'));d=mujoco.MjData(m);mujoco.mj_forward(m,d);bodies={}
    for g in range(m.ngeom):
        b=int(m.geom_bodyid[g]);name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g)
        if not b or m.geom_group[g]==3 or name in ['floor','v8_thorax_contact_band','elbow_joint']:continue
        if mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,b)=='thorax':continue
        solid,_=local_solid(m,d,g)
        if solid is not None:bodies.setdefault(b,[]).append(solid)
    rows=[]
    for b,solids in bodies.items():
        union=md.Manifold.batch_boolean(solids,md.OpType.Add);mesh=union.to_mesh();vertices=np.asarray(mesh.vert_properties)[:,:3]*.001;tri=vertices[np.asarray(mesh.tri_verts)]
        volumes=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2]))/6.;volume=volumes.sum();mass=density*volume;s=tri.sum(1)
        centroid=np.einsum('i,ij->j',volumes,s)/(4*volume)
        second=density*np.einsum('i,ijk->jk',volumes,(np.einsum('ij,ik->ijk',s,s)+np.einsum('imj,imk->ijk',tri,tri)))/20.
        inertia=np.trace(second)*np.eye(3)-second-mass*(np.dot(centroid,centroid)*np.eye(3)-np.outer(centroid,centroid))
        if mass<=0 or np.linalg.eigvalsh(inertia).min()<=0:raise RuntimeError('Invalid material-union mass properties')
        name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_BODY,b);body=root.find(f'.//body[@name="{name}"]');old=body.find('inertial')
        if old is not None:body.remove(old)
        ET.SubElement(body,'inertial',{'pos':v8.vec(centroid),'mass':str(mass),'fullinertia':v8.vec(np.r_[np.diag(inertia),inertia[0,1],inertia[0,2],inertia[1,2]])})
        rows.append({'body':name,'mass_kg':float(mass),'centroid_local_m':centroid.tolist(),'inertia_local_kg_m2':inertia.tolist()})
    (EXPORT/'prototype_mass_properties.json').write_text(json.dumps({'assumption':'Uniform 1100 kg/m3 printed material; union geometry; an estimate, not a weighed prototype; no payload','moving_part_mass_kg':sum(r['mass_kg'] for r in rows),'bodies':rows},indent=2))
