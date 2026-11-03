"""Single-source static wire paths for S1 interface design, not motion proof.

Bone centres and attachment seats are preserved. The shared deltoid cap has
three distinct candidate wire/crimp ports within its existing footprint.
No route interior is claimed as a fixed world guide in a moving model.
"""
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import frames
import shoulder_pads
import mesh_geom
import interference as I
import build_shoulder_tendon_segments as T
import search_shoulder_channel_paths as SP

DESIGN=frames.OUTPUT_DIR.parent/'shoulder_design'
DEST=DESIGN/'s1_port_path_revision'


def unit(v):
    v=np.asarray(v,float);n=np.linalg.norm(v)
    if n<1e-9:raise ValueError('Zero direction')
    return v/n


def trim_path(points,start,end):
    p=np.array(points,float);lengths=np.linalg.norm(np.diff(p,axis=0),axis=1)
    arc=np.r_[0,np.cumsum(lengths)]
    # np.sum and sequential cumsum can differ at the final few floating-point
    # bits. Accept only that arithmetic uncertainty, not a longer path.
    tol=2*np.finfo(float).eps*len(p)*max(1.,arc[-1])
    if -tol<=start<0:start=0.
    if arc[-1]<end<=arc[-1]+tol:end=float(arc[-1])
    if not 0<=start<end<=arc[-1]:raise ValueError('Invalid arclength interval')
    def at(s):
        i=min(np.searchsorted(arc,s,side='right')-1,len(p)-2)
        return p[i]+(p[i+1]-p[i])*(s-arc[i])/lengths[i]
    return np.vstack([at(start),p[(arc>start)&(arc<end)],at(end)])


def segment_distance(p,q,a,b):
    p,q,a,b=map(lambda x:np.asarray(x,float),(p,q,a,b));u=q-p;v=b-a;w=p-a
    aa=u@u;bb=u@v;cc=v@v;dd=u@w;ee=v@w
    candidates=[]
    for x in [p,q]:
        t=np.clip(((x-a)@v)/cc,0,1) if cc>1e-20 else 0
        candidates.append(np.linalg.norm(x-a-t*v))
    for x in [a,b]:
        s=np.clip(((x-p)@u)/aa,0,1) if aa>1e-20 else 0
        candidates.append(np.linalg.norm(x-p-s*u))
    det=aa*cc-bb*bb
    if det>1e-16:
        s=(bb*ee-cc*dd)/det;t=(aa*ee-bb*dd)/det
        if 0<=s<=1 and 0<=t<=1:candidates.append(np.linalg.norm(w+s*u-t*v))
    return float(min(candidates))


def path_distance(p,q):
    return min(segment_distance(a,b,c,d) for a,b in zip(p,p[1:]) for c,d in zip(q,q[1:]))


def outer_envelope(path,r=1.8):
    parts=[T.rod(a,b,r) for a,b in zip(path,path[1:]) if np.linalg.norm(b-a)>1e-6]
    parts += [I._manifold().Manifold.sphere(r,32).translate(p) for p in path]
    return I._manifold().Manifold.batch_boolean(parts,I._manifold().OpType.Add)


def main():
    DEST.mkdir(exist_ok=True)
    files=[DESIGN/'s1_modules.json',DESIGN/'s1_end_pads.json',
           DESIGN/'s1_integrated_revision/tendons_export.json']
    design=json.loads(files[0].read_text(encoding='utf8'))
    seats=json.loads(files[1].read_text(encoding='utf8'))['seats']
    tendons=json.loads(files[2].read_text(encoding='utf8'))['segments']
    seat_map={h:s for s in seats for h in s['heads']}
    guide_map={s['head']+'/'+s['end']:s for s in tendons if s['built']}
    bones={k:SP.load_solid(T.BONE_DIR/v) for k,v in SP.BONES.items()}
    tris=[]
    for v in SP.BONES.values():
        V,F=mesh_geom.load_obj(str(T.BONE_DIR/v));W,G=mesh_geom.weld(V,F);tris.append(W[G])
    cores={};labels={}
    for module in design['modules']:
        channels={r['head']:r for r in module['channels']}
        for row in module['corridors']:
            h=row['head'];p=np.array(row['path_source_mm']);total=np.linalg.norm(np.diff(p,axis=0),axis=1).sum()
            cut0=min(.35*total,max(12.,(channels[h].get('pad_to_channel_origin_mm') or 0)+3.))
            cut1=min(.35*total,max(12.,(channels[h].get('pad_to_channel_distal_mm') or 0)+3.))
            cores[h]=trim_path(p,cut0,total-cut1);labels[h]=module['module']
    shared=next(s for s in seats if len(s['heads'])==3)
    basis=shoulder_pads.pad_basis(shared['normal_source'])
    offsets=np.array([basis[:,0]*4.5*np.cos(a)+basis[:,1]*4.5*np.sin(a)
                      for a in np.radians([30,150,270])])
    # Pick the least-crossing orientation among the six assignments, without
    # changing the bone seat or recomputing an attachment from anatomy proximity.
    heads=[h.split('/')[0] for h in shared['heads']];centre=np.array(shared['centre_source_mm'])
    assignments=list(itertools.permutations(range(3)))
    chosen=min(assignments,key=lambda ids:sum(np.linalg.norm(cores[h][-1]-centre-offsets[j]) for h,j in zip(heads,ids)))
    shared_offsets={h:offsets[i] for h,i in zip(heads,chosen)}
    routes=[]
    report=dict(scope='static candidate unified wire-path specification; not manufacturing or actuation acceptance',
        source_sha256={str(p.relative_to(frames.PROJECT_ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        shared_port_radius_mm=4.5,wire_diameter_mm=.6,nominal_lumen_radius_mm=.85,
        rigid_guide_outer_radius_mm=1.8,soft_transition_outer_radius_mm=3.2,routes=routes,stage_pass=False)
    for h,core in cores.items():
        connectors={};ports={}
        for end,endpoint in [('origin',core[0]),('distal',core[-1])]:
            key=h+'/'+end;s=seat_map[key];n=unit(s['normal_source']);c=np.array(s['centre_source_mm'])
            offset=shared_offsets[h] if len(s['heads'])==3 else np.zeros(3)
            height=s['lift_mm']+s['pad_thickness_mm']+float(s.get('pad_boss_mm',3.))
            inner=c+offset+n*(height+.2);outer=inner+n*s['cap_thickness_mm']
            bend=inner-n*2
            old=guide_map.get(key)
            mid=np.array(old['path_source_mm'][1:]) if old else np.empty((0,3))
            initial=np.vstack([bend,mid,endpoint])
            clean=[initial[0]]+[b for a,b in zip(initial,initial[1:]) if np.linalg.norm(b-a)>1e-5]
            relaxed=T.relax_off_bone(clean,tris,3.5,step_mm=1.,iterations=50,lam=.2)
            path=np.vstack([inner,relaxed]);env=outer_envelope(path,r=3.2)
            hits={k:I.intersection_volume(env,b) for k,b in bones.items()}
            ports[end]=dict(bone_seat_source_mm=c.tolist(),offset_from_seat_mm=offset.tolist(),
                normal_source=n.tolist(),inner_face_source_mm=inner.tolist(),outer_face_source_mm=outer.tolist(),
                hand_pull_lead_source_mm=(outer+n*5).tolist() if end=='origin' else None)
            connectors[end]=dict(path_source_mm=path.tolist(),outer_envelope_bone_intersection_mm3=hits,
                bone_clear_at_nominal_geometry=all(v<=.001 for v in hits.values()))
            print(h,end,'connector bone envelope hits', {k:round(v,4) for k,v in hits.items()},flush=True)
        full=np.vstack([ports['origin']['hand_pull_lead_source_mm'],ports['origin']['outer_face_source_mm'],
            connectors['origin']['path_source_mm'],core[1:],
            np.array(connectors['distal']['path_source_mm'])[-2::-1],ports['distal']['outer_face_source_mm']])
        routes.append(dict(head=h,module=labels[h],core_path_source_mm=core.tolist(),
            ports=ports,connectors=connectors,full_wire_path_source_mm=full.tolist(),
            status='candidate_requires_clearance_and_material_generation'))
        (DEST/'routes_partial.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    deltoid=[r for r in routes if r['module']=='deltoid']
    report['deltoid_wire_pair_minimum_mm']=[dict(a=a['head'],b=b['head'],
        distance_mm=path_distance(np.array(a['full_wire_path_source_mm']),np.array(b['full_wire_path_source_mm'])))
        for a,b in itertools.combinations(deltoid,2)]
    report['missing']=['generated cap, pad, guide and TPU geometry from these routes',
        'neighbour-material clearance','wire bend limits','full joint constraints','moving shapes and measured properties']
    (DEST/'routes.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Static unified route specification saved; not released.',flush=True)


if __name__=='__main__':main()
