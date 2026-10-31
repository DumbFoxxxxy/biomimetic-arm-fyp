"""Regenerate neutral extensor geometry for the selected real attachment search.

This changes a manufactured muscle shape, not a loaded deformation result.
"""
from build_direct_muscle_yoke import *
import build_shaped_elbow_muscles as generator

def run():
    folder=OUT/'direct_muscle_interfaces/extensor_reanchor_extended'
    search=folder/'endpoint_search.json';result=json.loads(search.read_text());choice=result['selected'];assert choice is not None
    source=OUT/'elbow_muscles/FJ1477_profile.json';profile=json.loads(source.read_text());F=np.array(profile['basis_columns_display']);delta=np.array(choice['delta_display_mm']);local=F.T@delta;heights=np.array([s['height_mm'] for s in profile['sections']]);lo,hi=heights[[0,-1]];changes=[]
    for i,(h,s) in enumerate(zip(heights,profile['sections'])):
        t=np.clip(((h-lo)/(hi-lo)-.55)/.45,0,1);w=float(6*t**5-15*t**4+10*t**3)
        s['height_mm']=float(h+local[2]*w);s['boundary_local_mm']=(np.array(s['boundary_local_mm'])+local[:2]*w).tolist();s['pole_local_mm']=(np.array(s['pole_local_mm'])+local[:2]*w).tolist();changes.append({'section':i,'weight':w,'delta_display_mm':(delta*w).tolist()})
    updated=np.array([s['height_mm'] for s in profile['sections']]);assert np.all(np.diff(updated)>0)
    profile['body_height_range_mm']=[float(updated[0]),float(updated[-1])];profile['body_axial_span_mm']=float(updated[-1]-updated[0]);profile['engineering_revision']={'scope':'manufactured neutral shape and endpoint relocation, not compliant analysis','source_profile_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'endpoint_search_sha256':hashlib.sha256(search.read_bytes()).hexdigest(),'delta_display_mm':delta.tolist(),'proximal_boundary_unchanged':True,'bone_centres_changed':False,'sections':changes}
    dest=folder/'muscle';dest.mkdir(exist_ok=True);(dest/'FJ1477_profile.json').write_text(json.dumps(profile,indent=2),encoding='utf8');generator.DEST=dest;generator.run(('FJ1477',))

if __name__=='__main__':run()
