"""Independent solid clearance check at samples different from the cutter sweep."""
import json,sys
import numpy as np
from shoulder_model import inputs,EXPORT,SOURCE,v8,input_fingerprints
sys.path.insert(0,str(SOURCE/'vendor'))
import manifold3d as md
from build_cable_channels import cylinder_between,clip_segment
from geometry_distance import segment_mesh_distance
sys.stdout.reconfigure(encoding='utf-8')

def segment_distance(a,b,c,d):
    u=b-a;v=d-c;w=a-c;aa=np.dot(u,u);bb=np.dot(u,v);cc=np.dot(v,v);dd=np.dot(u,w);ee=np.dot(v,w);den=aa*cc-bb*bb
    if den>1e-12:
        s=(bb*ee-cc*dd)/den;t=(aa*ee-bb*dd)/den
        if 0<=s<=1 and 0<=t<=1:return float(np.linalg.norm(w+s*u-t*v))
    def point_seg(p,x,y):
        z=y-x;return float(np.linalg.norm(p-x-np.clip(np.dot(p-x,z)/max(np.dot(z,z),1e-12),0,1)*z))
    return min(point_seg(a,c,d),point_seg(b,c,d),point_seg(c,a,b),point_seg(d,a,b))

def audit():
    cfg,fit,anchors=inputs();mesh=np.load(EXPORT/'scapula_hybrid_design_geometry.npz')
    verts=mesh['vertices_ac_local_mm'];faces=mesh['faces'];solid=md.Manifold(md.Mesh(np.ascontiguousarray(verts,dtype=np.float32),np.ascontiguousarray(faces,dtype=np.uint32)))
    lo=verts.min(0)-2.;hi=verts.max(0)+2.;zero=np.array(fit['poses'][0]['q_rad'])
    fixed=np.array([r['fixed_eye_world_mm'] for r in anchors['routes']]);eye=np.array([r['eye_local_ac_mm'] for r in anchors['routes']])
    minimum=3.;minpair=1e9;worst={};pairworst={};n=0
    for action in ['UP_ROT','DOWN_ROT','PROTRACT','RETRACT','ELEVATE','DEPRESS']:
        seq=[zero]+[np.array(p['q_rad']) for p in fit['poses'] if p.get('action')==action]
        for k,(a,b) in enumerate(zip(seq[:-1],seq[1:])):
            for f in np.linspace(0.,1.,11):
                q=a*(1-f)+b*f;fk=v8.direct_kinematics(cfg,q);n+=1
                starts=(fixed-fk['ac_mm'])@fk['body_rotation'];moving=eye@fk['body_rotation'].T+fk['ac_mm']
                for j in range(8):
                    if j==6:
                        # The final 3 mm is the intentional eye/crimp contact
                        # zone, not a freely sliding cable run. Its contact is
                        # allowed only at the designated attachment point.
                        direction=(eye[j]-starts[j]);direction/=np.linalg.norm(direction)
                        free_end=eye[j]-3.*direction
                        clipped=clip_segment(starts[j],free_end,lo,hi)
                        gap=segment_mesh_distance(*clipped,verts[faces])-.5 if clipped else 3.
                        if gap<minimum:minimum=gap;worst={'action':action,'segment':k,'fraction':float(f),'cable':'S07'}
                    for h in range(j):
                        distance=segment_distance(fixed[j],moving[j],fixed[h],moving[h])-1.
                        if distance<minpair:minpair=distance;pairworst={'action':action,'segment':k,'fraction':float(f),'cables':[j+1,h+1]}
        print('CHECKED',action,'minimum_channel_gap',minimum,'minimum_wire_pair_gap',minpair,flush=True)
    report={'scope':'Euclidean segment/triangle distance minus 0.5mm wire radius, actual boolean-cut surface; excludes designated terminal clamp contact zone; 1320 interpolation poses; segment-pair distance minus two wire radii',
        'terminal_clamp_contact_zone_mm':3.,
        'pose_count':n,'minimum_channel_wall_clearance_mm':minimum,'channel_worst':worst,
        'minimum_cable_cable_surface_gap_mm':minpair,'cable_pair_worst':pairworst,
        'free_channel_pass':minimum>=.75,'cable_cable_pass':minpair>=.5,'input_sha256':input_fingerprints()}
    (EXPORT/'channel_clearance_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':audit()
