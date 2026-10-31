"""Anatomy-derived sectional design coordinates, not inferred muscle fibers."""
from repair_inputs import *
from scipy.interpolate import make_interp_spline
from build_SC_connection import unit
DEST=OUT/'elbow_muscles'

def area(poly):
    q=np.roll(poly,-1,axis=0)
    return float(np.sum(poly[:,0]*q[:,1]-poly[:,1]*q[:,0])/2)

def inside(points,poly):
    result=np.zeros(len(points),dtype=bool)
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        crossing=(a[1]>points[:,1])!=(b[1]>points[:,1])
        result^=crossing&(points[:,0]<(b[0]-a[0])*(points[:,1]-a[1])/(b[1]-a[1]+1e-30)+a[0])
    return result

def boundary_distance(points,poly):
    distance=np.full(len(points),np.inf)
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        d=b-a;t=np.clip((points-a)@d/(d@d+1e-30),0,1)
        distance=np.minimum(distance,np.linalg.norm(points-a-t[:,None]*d,axis=1))
    return distance

def section(v,f,height):
    tri=v[f];segments=[]
    for t in tri:
        hits=[]
        for a,b in zip(t,np.roll(t,-1,axis=0)):
            if (a[2]-height)*(b[2]-height)<0:
                hits.append((a+(b-a)*(height-a[2])/(b[2]-a[2]))[:2])
        if len(hits)==2:segments.append(hits)
    nodes={};positions=[];edges=[]
    for segment in segments:
        pair=[]
        for p in segment:
            key=tuple(np.round(p,6))
            if key not in nodes:nodes[key]=len(positions);positions.append(p)
            pair.append(nodes[key])
        if pair[0]!=pair[1]:edges.append(pair)
    graph={i:[] for i in range(len(positions))}
    for a,b in edges:graph[a].append(b);graph[b].append(a)
    if any(len(a)!=2 for a in graph.values()):raise RuntimeError('Non-simple section graph')
    remaining=set(graph);loops=[]
    while remaining:
        start=next(iter(remaining));current=start;previous=None;loop=[]
        while True:
            loop.append(current);remaining.discard(current);choices=graph[current];nxt=choices[0] if choices[0]!=previous else choices[1]
            previous,current=current,nxt
            if current==start:break
        poly=np.array(positions)[loop]
        if area(poly)<0:poly=poly[::-1]
        loops.append(poly)
    return sorted(loops,key=lambda p:abs(area(p)),reverse=True)

def resample(poly,n=32):
    # Common cross-sectional +X extremum defines a reproducible seam.
    poly=np.roll(poly,-int(np.argmax(poly[:,0])),axis=0);closed=np.vstack([poly,poly[0]])
    lengths=np.r_[0,np.cumsum(np.linalg.norm(np.diff(closed,axis=0),axis=1))]
    q=np.linspace(0,lengths[-1],n,endpoint=False)
    return np.column_stack([np.interp(q,lengths,closed[:,i]) for i in [0,1]])

def kernel(poly):
    lo=poly.min(0)-1;hi=poly.max(0)+1
    result=np.array([[lo[0],lo[1]],[hi[0],lo[1]],[hi[0],hi[1]],[lo[0],hi[1]]])
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        new=[];d=b-a
        if len(result)<3:return np.empty((0,2))
        signed=lambda p:d[0]*(p[1]-a[1])-d[1]*(p[0]-a[0])
        for p,q in zip(result,np.roll(result,-1,axis=0)):
            sp,sq=signed(p),signed(q)
            if sp>=-1e-8:new.append(p)
            if (sp>0)!=(sq>0):new.append(p+(q-p)*sp/(sp-sq))
        result=np.array(new)
    return result

def pole(poly):
    k=kernel(poly)
    if len(k)>=3:
        point=k.mean(0)
        if boundary_distance(point[None],poly)[0]>1:return point,True
    x=np.linspace(*[poly[:,0].min(),poly[:,0].max()],45);y=np.linspace(*[poly[:,1].min(),poly[:,1].max()],45)
    grid=np.array(list(itertools.product(x,y)));grid=grid[inside(grid,poly)]
    if not len(grid):raise RuntimeError('No interior section pole')
    dist=boundary_distance(grid,poly);return grid[np.argmax(dist)],False

def route_engineering_poles(sections):
    """Shortest-step interior path with fixed end poles; NOT anatomical fibers."""
    candidates=[];clearances=[]
    for i,s in enumerate(sections):
        poly=np.array(s['boundary_local_mm']);original=np.array(s['pole_local_mm'])
        if i in [0,len(sections)-1]:points=original[None]
        else:
            x=np.linspace(poly[:,0].min(),poly[:,0].max(),35);y=np.linspace(poly[:,1].min(),poly[:,1].max(),35);points=np.array(list(itertools.product(x,y)));points=points[inside(points,poly)]
            distance=boundary_distance(points,poly);points=points[distance>=1.]
            if not len(points):points=original[None]
            if len(points)>220:points=points[np.unique(np.linspace(0,len(points)-1,220).astype(int))]
            points=np.vstack([points,original])
        candidates.append(points);clearances.append(boundary_distance(points,poly))
    cost=np.zeros(1);parents=[]
    for i in range(1,len(candidates)):
        delta=candidates[i][:,None,:]-candidates[i-1][None,:,:]
        values=(delta*delta).sum(2)+cost[None,:]
        parent=values.argmin(1);cost=values[np.arange(len(parent)),parent]+.5/(clearances[i]+.1)**2;parents.append(parent)
    index=int(cost.argmin());result=[candidates[-1][index]]
    for i in range(len(sections)-2,-1,-1):index=int(parents[i][index]);result.append(candidates[i][index])
    return np.array(list(reversed(result)))

def run():
    DEST.mkdir(exist_ok=True);definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'))
    E=np.array(definition['centres_source_mm']['elbow'])-SHIFT;G=np.array(definition['centres_source_mm']['GH'])-SHIFT
    z=unit(E-G);x=unit(np.cross(z,[0,1,0]));y=np.cross(z,x);F=np.stack([x,y,z],1);rows=[]
    for fid in ['FJ1486','FJ1477']:
        path=OUT/'repaired'/(fid+'.obj');v,f=read_obj(str(path));local=(v-SHIFT-E)@F
        low,high=local[:,2].min(),local[:,2].max();heights=np.linspace(low+.10*(high-low),high-.10*(high-low),17);sections=[]
        for height in heights:
            loops=section(local,f,height);poly=loops[0];centre,star=pole(poly);sample=resample(poly)
            sections.append({'height_mm':float(height),'boundary_local_mm':sample.tolist(),'pole_local_mm':centre.tolist(),'pole_boundary_gap_mm':float(boundary_distance(centre[None],poly)[0]),'pole_in_visibility_kernel':star,'discarded_secondary_loop_areas_mm2':[abs(area(p)) for p in loops[1:]],'source_section_area_mm2':abs(area(poly))})
        row={'id':fid,'source':str(path.relative_to(ROOT)),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'coordinate_origin_display_mm':E.tolist(),'basis_columns_display':F.tolist(),'source_height_range_mm':[float(low),float(high)],'body_height_range_mm':[float(heights[0]),float(heights[-1])],'body_axial_span_mm':float(heights[-1]-heights[0]),'sections':sections,'engineering_axis':'humeral proximal-to-distal axis; transverse frame fixed in source coordinates, NOT PCA/fiber inference','surface_design':'sectionwise periodic cubic B-spline, cubic interpolation along axial parameter; weights=1, a nonrational NURBS subclass','fiber_direction_verified':False,'anatomical_end_cropping_fraction':.1,'frozen':False}
        (DEST/(fid+'_profile.json')).write_text(json.dumps(row,indent=2),encoding='utf8');rows.append(row)
        print(fid,'span',round(row['body_axial_span_mm'],2),'min pole gap',round(min(s['pole_boundary_gap_mm'] for s in sections),2),'star sections',sum(s['pole_in_visibility_kernel'] for s in sections),'of',len(sections),flush=True)
    (DEST/'profiles.json').write_text(json.dumps({'scope':'anatomical section design inputs, not final muscle or deformation acceptance','muscles':[r['id'] for r in rows]},indent=2),encoding='utf8')

if __name__=='__main__':run()
