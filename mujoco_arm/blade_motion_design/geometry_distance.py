"""Independent Euclidean segment-to-triangle distance (all units consistent)."""
import numpy as np

def point_segments_sq(p,a,b):
    d=b-a;t=np.sum((p-a)*d,axis=-1)/np.maximum(np.sum(d*d,axis=-1),1e-24)
    delta=p-a-np.clip(t,0,1)[...,None]*d
    return np.sum(delta*delta,axis=-1)

def segment_mesh_distance(a,b,triangles):
    triangles=np.asarray(triangles,dtype=np.float64);a=np.asarray(a,dtype=np.float64);b=np.asarray(b,dtype=np.float64)
    A=triangles[:,0];B=triangles[:,1];C=triangles[:,2];ab=B-A;ac=C-A
    direction=b-a;h=np.cross(direction,ac);det=np.sum(ab*h,axis=1);ok=np.abs(det)>1e-12
    inv=np.where(ok,1/np.where(ok,det,1),0.);s=a-A;uu=inv*np.sum(s*h,axis=1);q=np.cross(s,ab)
    vv=inv*np.sum(direction*q,axis=1);tt=inv*np.sum(ac*q,axis=1)
    if np.any(ok&(uu>=0)&(vv>=0)&(uu+vv<=1)&(tt>=0)&(tt<=1)):return 0.
    # Endpoints whose orthogonal projection is inside a triangle.
    normal=np.cross(ab,ac);n2=np.sum(normal*normal,axis=1)
    d00=np.sum(ab*ab,axis=1);d01=np.sum(ab*ac,axis=1);d11=np.sum(ac*ac,axis=1);den=d00*d11-d01*d01
    candidates=[]
    for p in [a,b]:
        signed=np.sum((p-A)*normal,axis=1);proj=p-normal*(signed/np.maximum(n2,1e-24))[:,None];v2=proj-A
        d20=np.sum(v2*ab,axis=1);d21=np.sum(v2*ac,axis=1)
        u=(d11*d20-d01*d21)/np.maximum(den,1e-24);v=(d00*d21-d01*d20)/np.maximum(den,1e-24)
        inside=(den>1e-16)&(u>=0)&(v>=0)&(u+v<=1)
        if np.any(inside):candidates.append(float(np.min(signed[inside]**2/n2[inside])))
    # Interior segment/edge minima, plus every boundary case.
    for p,z in [(A,B),(B,C),(C,A)]:
        u=direction;v=z-p;w=a-p;aa=np.dot(u,u);bb=v@u;cc=np.sum(v*v,axis=1);dd=w@u;ee=np.sum(v*w,axis=1);de=aa*cc-bb*bb
        safe=np.where(np.abs(de)>1e-20,de,1.)
        s=(bb*ee-cc*dd)/safe;t=(aa*ee-bb*dd)/safe
        valid=(np.abs(de)>1e-20)&(s>=0)&(s<=1)&(t>=0)&(t<=1)
        if np.any(valid):
            delta=w[valid]+s[valid,None]*u-t[valid,None]*v[valid]
            candidates.append(float(np.min(np.sum(delta*delta,axis=1))))
        candidates.extend([float(np.min(point_segments_sq(a,p,z))),float(np.min(point_segments_sq(b,p,z))),
                           float(np.min(point_segments_sq(p,a,b))),float(np.min(point_segments_sq(z,a,b)))])
    return float(np.sqrt(max(0.,min(candidates))))
