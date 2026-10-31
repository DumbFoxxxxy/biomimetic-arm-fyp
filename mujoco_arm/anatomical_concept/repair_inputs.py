"""Close only simple triangular holes, preserving original vertex coordinates.
Repair is explicit computational preprocessing, not a claim of anatomy recovery.
"""
from build_reference import *
sys.path.insert(0,str(ROOT/'tmp/shoulder_motion_design_20261029/vendor'))
import manifold3d as mf

def run():
    supports='--supports' in sys.argv
    report=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));dest=OUT/('repaired_supports' if supports else 'repaired');dest.mkdir(exist_ok=True);rows=[]
    for record in report['models']:
        if record['kind']!=('support_reference' if supports else 'muscle'):continue
        v,f=read_obj(str(ROOT/record['source']));v,f=weld_stl(v[f]);edge_faces={};small_open=[]
        parent=list(range(len(v)))
        def find(a):
            while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
            return a
        for face in f:
            for b in face[1:]:parent[find(int(b))]=find(int(face[0]))
        components={}
        for i,face in enumerate(f):components.setdefault(find(int(face[0])),[]).append(i)
        remove_faces=set()
        for indices in components.values():
            if not 1<len(indices)<=4:continue
            pts=v[np.unique(f[indices])];span=float(np.linalg.norm(pts.max(0)-pts.min(0)));edges=Counter(tuple(sorted((int(a),int(b)))) for face in f[indices] for a,b in zip(face,np.roll(face,-1)))
            if span<1 and any(n==1 for n in edges.values()):
                small_open.append({'original_face_indices':indices,'max_span_mm':span,'reason':'detached open surface component, <=4 faces and <1mm span; not a closed tissue volume'});remove_faces.update(indices)
        if remove_faces:f=f[[i for i in range(len(f)) if i not in remove_faces]]
        for face in f:
            for a,b in zip(face,np.roll(face,-1)):edge_faces.setdefault(tuple(sorted((int(a),int(b)))),[]).append((int(a),int(b)))
        directed=[p[0] for p in edge_faces.values() if len(p)==1];lookup={a:b for a,b in directed};remaining=set(lookup);patches=[];failed=[]
        while remaining:
            start=next(iter(remaining));cycle=[start];cur=lookup[start]
            while cur!=start and cur not in cycle and len(cycle)<=len(lookup):cycle.append(cur);cur=lookup.get(cur,start)
            remaining-=set(cycle)
            if len(cycle)!=3 or cur!=start:failed.append(cycle);continue
            pts=v[cycle];area=float(np.linalg.norm(np.cross(pts[1]-pts[0],pts[2]-pts[0]))/2);span=max(float(np.linalg.norm(a-b)) for a,b in itertools.combinations(pts,2));patches.append({'vertices':cycle,'area_mm2':area,'max_span_mm':span})
        # Some triangular boundaries are detached single-face islands, not holes.
        # Capping these creates zero-thickness double faces and fails OBJ roundtrip.
        keys={tuple(sorted(face)) for face in f};isolated=[p for p in patches if tuple(sorted(p['vertices'])) in keys];patches=[p for p in patches if tuple(sorted(p['vertices'])) not in keys]
        removed={tuple(sorted(p['vertices'])) for p in isolated};faces=[face.tolist() for face in f if tuple(sorted(face)) not in removed]+[list(reversed(p['vertices'])) for p in patches]
        solid=mf.Manifold(mf.Mesh64(v.astype(np.float64),np.array(faces,dtype=np.uint64)))
        status=str(solid.status());volume=float(solid.volume()) if solid.status()==mf.Error.NoError else None
        row={'id':record['id'],'name':record['name'],'added_triangles':len(patches),'patches':patches,'removed_small_open_surface_components':small_open,'removed_zero_thickness_triangle_islands':isolated,'unresolved_boundary_components':failed,'manifold_status':status,'volume_mm3':volume,'vertices_moved':False,'clinical_attachment_inference':False}
        if solid.status()==mf.Error.NoError and not failed:
            if volume<0:faces=[face[::-1] for face in faces];row['global_winding_reversed']=True
            # OBJ retains full source precision and shares the raw bone frame.
            text=['# Computationally cleaned copy; see repair_report.json for all changes.','# BodyParts3D, (c) The Database Center for Life Science licensed under CC Attribution-Share Alike 2.1 Japan.','# http://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html']
            text+=['v '+' '.join(map(str,p)) for p in v];text+=['f '+' '.join(str(int(i)+1) for i in face) for face in faces]
            (dest/(record['id']+'.obj')).write_text('\n'.join(text)+'\n',encoding='utf8')
        rows.append(row)
    data={'scope':'detached zero-thickness face islands removed; actual triangular holes capped only when no existing face occupies the loop; no vertex motion or clinical reconstruction','source_files_untouched':True,'models' if supports else 'muscles':rows};(OUT/('support_repair_report.json' if supports else 'repair_report.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf8');print([(r['id'],len(r['removed_zero_thickness_triangle_islands']),r['manifold_status']) for r in rows],flush=True)

if __name__=='__main__':run()
