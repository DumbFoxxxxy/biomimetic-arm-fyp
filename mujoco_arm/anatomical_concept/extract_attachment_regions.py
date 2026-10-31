"""Atlas-constrained bone proximity patches; never automatic clinical landmarks.
Exact point-to-triangle distances at muscle FACE CENTROIDS. Areas are sampled
proximity estimates, not verified physical attachment areas or drill locations.
"""
from build_reference import *
from attachment_atlas import atlas,BONES
from scipy.spatial import cKDTree

def closest(p,t):
    a,b,c=t[:,0],t[:,1],t[:,2];u=b-a;v=c-a;n=np.cross(u,v);nn=np.sum(n*n,axis=1);s=np.sum((p-a)*n,axis=1);proj=p-n*(s/np.maximum(nn,1e-30))[:,None];w=proj-a;uu=np.sum(u*u,axis=1);uv=np.sum(u*v,axis=1);vv=np.sum(v*v,axis=1);wu=np.sum(w*u,axis=1);wv=np.sum(w*v,axis=1);den=uu*vv-uv*uv;alpha=(vv*wu-uv*wv)/np.maximum(den,1e-30);beta=(uu*wv-uv*wu)/np.maximum(den,1e-30);inside=(den>1e-20)&(alpha>=0)&(beta>=0)&(alpha+beta<=1)
    candidates=[proj];errors=[np.where(inside,np.sum((proj-p)**2,axis=1),np.inf)]
    for x,y in [(a,b),(b,c),(c,a)]:
        d=y-x;r=np.clip(np.sum((p-x)*d,axis=1)/np.maximum(np.sum(d*d,axis=1),1e-30),0,1);q=x+r[:,None]*d;candidates.append(q);errors.append(np.sum((q-p)**2,axis=1))
    ei,ti=np.unravel_index(np.argmin(np.array(errors)),(4,len(t)));return candidates[ei][ti],float(np.sqrt(errors[ei][ti])),int(ti)

def projector(v,f):
    tri=v[f];c=tri.mean(1);radius=float(np.max(np.linalg.norm(tri-c[:,None],axis=2)));tree=cKDTree(c)
    def project(points):
        _,first=tree.query(points);positions=[];distances=[];indices=[]
        for p,i in zip(points,first):
            _,upper,_=closest(p,tri[[i]]);potential=tree.query_ball_point(p,upper+radius+1e-7);q,d,j=closest(p,tri[potential]);positions.append(q);distances.append(d);indices.append(potential[j])
        return np.array(positions),np.array(distances),np.array(indices)
    return project

def run():
    inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));known=atlas();dest=OUT/'attachment_regions';dest.mkdir(exist_ok=True);bones={};sources={r['id']:r for r in inv['models']};reference_ids={name:[fid] for name,fid in BONES.items()}
    for name,fid in BONES.items():
        v,f=read_obj(str(ROOT/sources[fid]['source']));v,f=weld_stl(v[f]);bones[name]=projector(v,f)
    support_file=ROOT/'骨骼/支撑参考/sources.json'
    if support_file.exists():
        support=json.loads(support_file.read_text(encoding='utf8'));groups={}
        for rec in support['models']:
            for category in rec['categories']:groups.setdefault(category,[]).append(rec)
        for category,refs in groups.items():
            vertices=[];faces=[];offset=0
            for rec in refs:
                v,f=read_obj(str(ROOT/rec['workspace_source']));v,f=weld_stl(v[f]);vertices.append(v);faces.append(f+offset);offset+=len(v)
            bones[category]=projector(np.concatenate(vertices),np.concatenate(faces));reference_ids[category]=[r['id'] for r in refs]
    rows=[];tree=ET.parse(OUT/'reference.xml');root=tree.getroot();world=root.find('worldbody')
    for rec in inv['models']:
        if rec['kind']!='muscle':continue
        fid=rec['id']
        if fid not in known:raise RuntimeError('Unmapped muscle '+fid)
        spec=known[fid];v,f=read_obj(str(OUT/'repaired'/(fid+'.obj')));tri=v[f];centres=tri.mean(1);areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2;regions=[];projection_cache={}
        for role in ['origin','insertion']:
            for structure in spec[role+'_structures']:
                if structure not in bones:
                    regions.append({'role':role,'structure':structure,'status':'support geometry absent from current five bones; must design frame interface or add anatomical reference','coordinate_fabricated':False});continue
                if structure not in projection_cache:projection_cache[structure]=bones[structure](centres)
                points,distance,indices=projection_cache[structure];selected=np.where((distance<=2)&(areas>1e-8))[0];broad5=np.where((distance<=5)&(areas>1e-8))[0];key=fid+'_'+role+'_'+structure
                np.savez_compressed(dest/(key+'.npz'),muscle_face_indices=selected,muscle_centres_source_mm=centres[selected],bone_points_source_mm=points[selected],bone_triangle_indices=indices[selected],distances_mm=distance[selected],muscle_face_areas_mm2=areas[selected])
                regions.append({'role':role,'structure':structure,'target_reference_ids':reference_ids[structure],'named_region':spec[role+'_named_region'],'status':'unverified atlas-constrained proximity patch' if len(selected) else 'no <=2mm proximity patch; needs landmark review','selected_face_count_2mm':len(selected),'selected_face_count_5mm':len(broad5),'sampled_muscle_surface_area_within_2mm_mm2':float(areas[selected].sum()),'minimum_centroid_gap_mm':float(distance.min()),'data':key+'.npz','semantic_landmark_localized':False,'centroid_is_not_drill_point':True})
                if len(selected):
                    draw=selected[np.unique(np.linspace(0,len(selected)-1,min(20,len(selected))).astype(int))]
                    for j,idx in enumerate(draw):ET.SubElement(world,'site',name=key+'_'+str(j),pos=' '.join(str(c/1000) for c in points[idx]-SHIFT),size='.0015',rgba='0 .9 .2 1' if role=='origin' else '1 .65 0 1',group='3')
        rows.append({'id':fid,'name':rec['name'],'anatomical_structure_constraints':spec,'regions':regions,'fiber_direction_verified':False,'frozen':False});print(fid,'patches',sum(r.get('selected_face_count_2mm',0)>0 for r in regions),'missing supports',sum('coordinate_fabricated' in r for r in regions),flush=True)
    ET.indent(root);root.find('compiler').set('meshdir','../assets');tree.write(dest/'reference.xml',encoding='utf8',xml_declaration=True)
    report={'scope':'27 anatomical muscle structure mappings and surface-distance probes on upper-limb and imported support references; no origin/insertion landmark or cable-anchor verification','method':'atlas determines allowed structures first; proximity never chooses anatomical identity. Exact nearest triangle is bounded using centroid distance plus maximum triangle radius. Patch threshold2mm, additional5mm count is diagnostic only.','muscles':rows,'missing_support_structure_types':sorted({r['structure'] for row in rows for r in row['regions'] if 'coordinate_fabricated' in r})};(dest/'atlas_and_regions.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    text=['# 肌肉附着关系与几何候选区域','','以下结构对应关系按大学解剖图谱整理。几何区域只是在这些指定骨骼上做表面距离筛查，','不等同于已确认的解剖附着区、钢丝锚点或可钻孔坐标。PCA方向没有用于推断肌纤维。','','|ID|肌肉|起始结构|止端结构|参考来源|','|---|---|---|---|---|']
    for r in rows:
        spec=r['anatomical_structure_constraints'];text.append(f"|{r['id']}|{r['name']}|{spec['origin_named_region']}|{spec['insertion_named_region']}|[atlas]({spec['source_url']})|")
    text+=['','## 尚缺的支撑定义','','当前仅有五块上肢骨骼。以下结构并没有被替换成虚构固定点：',', '.join(report['missing_support_structure_types']),'','后续应为这些区域建立真实参考或明确的工程机架连接。','斜方肌分部的整体附着范围仅用于候选筛查，不要求每个分部同时连接锁骨与肩胛骨。','原点和止端的骨骼对应关系确定后，还需在模型上定位具体骨性标志，并设计分布式锚定。','','NPZ保留所有被筛选的肌肉面、骨面投影点和距离。标记为绿色起始区域、橙色止端区域。','显示20个抽样点不意味着只有20个锚点，也不意味着所有候选点均可制造。'];(dest/'ATTACHMENT_REVIEW.md').write_text('\n'.join(text)+'\n',encoding='utf8')

if __name__=='__main__':run()
