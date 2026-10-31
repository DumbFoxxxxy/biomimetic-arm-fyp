"""Localize the first full-step adapter/scapula interference before rerouting."""
from build_elbow_muscle_mounts import *
from scipy.spatial.transform import Rotation

def run():
    data=json.loads((DEST/'assembly.json').read_text());paths=json.loads((DEST/'motion_paths.json').read_text());pose=next((p for p in paths['GH_axial_positive'] if any(set(h['pair'])=={'FJ3384','FJ1477_bone_adapter'} for h in p['hits'])),None)
    if not pose:print('No adapter limit on prescribed positive axial path');return
    definition=json.loads((OUT/'joint_model/definition.json').read_text());G=np.array(definition['centres_source_mm']['GH'])-SHIFT;R=Rotation.from_rotvec(np.radians(pose['controls_deg'][6:9])).as_matrix();adapter=load(DEST/'FJ1477_bone_adapter.obj')[0];scaprec=next(r for r in data['instances'] if r['key']=='FJ3384');scap=load(ROOT/scaprec['file'])[0]
    intersection=adapter.transform(np.c_[R,G-R@G])^scap;mesh=intersection.to_mesh64();points=np.asarray(mesh.vert_properties)[:,:3];point=(points.mean(0)-G)@R+G
    p=next(p for p in data['placements'] if p['id']=='FJ1477');segments=[]
    for i,path in enumerate(p['brace_paths_display_mm']):
        for j,(a,b) in enumerate(zip(np.array(path),np.array(path)[1:])):
            d=b-a;t=np.clip((point-a)@d/(d@d),0,1);segments.append({'branch':i,'segment':j,'centreline_distance_mm':float(np.linalg.norm(point-a-t*d))})
    segments.sort(key=lambda r:r['centreline_distance_mm']);report={'scope':'first full-degree colliding sample, not refined boundary','controls_deg':pose['controls_deg'],'contact_centroid_in_neutral_display_mm':point.tolist(),'overlap_mm3':float(intersection.volume()),'ranked_branch_segments':segments}
    (DEST/'adapter_limit_localization.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)

if __name__=='__main__':run()
