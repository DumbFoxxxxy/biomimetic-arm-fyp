"""All-pair rigid checks and sampled elbow paths for real M4 lever solids."""
from build_elbow_side_lever import *
from scipy.spatial.transform import Rotation
import itertools

def run():
    assembly=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'));baseline={};owners={};soft={}
    for r in assembly['instances']:
        s=load(ROOT/r['file'])[0].translate(r['translation_display_mm'])
        if r['owner'] in ['soft_unvalidated','distal_unvalidated']:soft[r['key']]=s;continue
        if r['key']=='FJ3391':continue
        baseline[r['key']]=s;owners[r['key']]=r['owner']
    support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in support['models']:
        if r['kind']=='bone_reference' and r['id'] not in baseline:
            baseline[r['id']]=load(OUT/'repaired_supports'/(r['id']+'.obj'),True)[0];owners[r['id']]='world'
    fixed={k:s for k,s in baseline.items() if owners[k]!='elbow'};existing_moving={k:s for k,s in baseline.items() if owners[k]=='elbow'}
    def hits(a,sa,b,sb):
        aa=np.array(sa.bounding_box());bb=np.array(sb.bounding_box());delta=np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])
        if np.any(delta>=0):return None
        volume=max(0,float((sa^sb).volume()))
        return {'pair':[a,b],'volume_mm3':volume} if volume>.001 else None
    results=[]
    variants=json.loads((DEST/'variants.json').read_text(encoding='utf8'))['variants']
    for variant in variants:
        folder=DEST/variant;design=json.loads((folder/'design.json').read_text(encoding='utf8'));parts={k:load(folder/(k+'.obj'))[0] for k in design['parts']};moving={**existing_moving,**parts};neutral=[];soft_flags=[]
        for (a,sa),(b,sb) in itertools.combinations(moving.items(),2):
            hit=hits(a,sa,b,sb)
            if hit:neutral.append(hit)
        for a,sa in parts.items():
            for b,sb in fixed.items():
                hit=hits(a,sa,b,sb)
                if hit:neutral.append(hit)
            for b,sb in soft.items():
                hit=hits(a,sa,b,sb)
                if hit:soft_flags.append(hit)
        path=[];back=[];E=np.array(design['joint_centre_display_mm']);axis=np.array(design['joint_axis'])
        def pose(angle):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();transform=np.c_[R,E-R@E];contacts=[];gap=3.;closest=[]
            for a,local in moving.items():
                moved=local.transform(transform);aa=np.array(moved.bounding_box())
                for b,sb in fixed.items():
                    bb=np.array(sb.bounding_box());delta=np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])
                    if np.linalg.norm(np.maximum(delta,0))>=3:continue
                    hit=hits(a,moved,b,sb)
                    if hit:contacts.append(hit)
                    value=float(moved.min_gap(sb,3))
                    if value<gap:gap=value;closest=[a,b]
            return {'angle_deg':float(angle),'hits':contacts,'moving_to_fixed_gap_mm_capped_3':gap,'closest_pair':closest}
        safe=0.;boundary=None
        if not neutral:
            path=[pose(angle) for angle in range(151)]
            first=next((i for i,p in enumerate(path) if p['hits']),None)
            if first is None:safe=150.
            elif first>0:
                low,high=first-1.,float(first)
                while high-low>.1:
                    mid=(low+high)/2
                    if pose(mid)['hits']:high=mid
                    else:low=mid
                safe=low;boundary=pose(high)
            else:boundary=path[0]
            back=[pose(angle) for angle in range(150,-1,-1)]
        report={'variant':variant,'manufactured_new_parts':len(parts),'new_parts_single_component':all(v==1 for v in design['parts'].values()),'moving_parts':list(moving),'fixed_count_including_bony_context':len(fixed),'neutral_rigid_interferences':neutral,'neutral_unvalidated_TPU_or_distal_contact_flags':soft_flags,'angular_step_deg':1,'boundary_threshold_deg':.1,'continuous_rigid_elbow_prefix_deg':safe,'first_motion_limit':boundary,'forward':path,'return':back,'rigid_elbow_path_clear':not neutral and safe==150,'not_included':['soft deformation or free distal end motion','tendon route and drive','shoulder13path audit','threads/preload/strength/printing'],'source_mesh_sha256':{k:hashlib.sha256((folder/(k+'.obj')).read_bytes()).hexdigest() for k in parts},'frozen':False}
        (folder/'geometry_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');results.append({k:report[k] for k in ['variant','neutral_rigid_interferences','continuous_rigid_elbow_prefix_deg','first_motion_limit','rigid_elbow_path_clear']})
        print(variant,'neutral',len(neutral),'softflags',len(soft_flags),'prefix',safe,'first',neutral[:4] if neutral else boundary,flush=True)
    (DEST/'comparison.json').write_text(json.dumps({'scope':'rigid elbow-side lever comparison only; no selected actuator/transmission','variants':results},indent=2),encoding='utf8')
if __name__=='__main__':run()
