"""Real anatomy-derived loft modification for the direct-attachment candidate.

Preserves proximal boundary and bones; smooth displacement only in distal45%.
This is a new neutral manufacturing shape, not a simulated loaded deformation.
"""
from build_forearm_receivers import ROOT,OUT,np
import json,hashlib
DEST=OUT/'direct_muscle_interfaces/shifted_muscle'

def run():
    DEST.mkdir(exist_ok=True);search=json.loads((OUT/'direct_muscle_interfaces/flexor_attachment_search.json').read_text());chosen=search['selected'];assert chosen is not None;delta=np.array(chosen['delta_display_mm']);source=OUT/'elbow_muscles/FJ1486_profile.json';profile=json.loads(source.read_text());F=np.array(profile['basis_columns_display']);local=F.T@delta;heights=np.array([r['height_mm'] for r in profile['sections']]);lo,hi=heights[[0,-1]];changes=[]
    for i,(h,s) in enumerate(zip(heights,profile['sections'])):
        t=np.clip(((h-lo)/(hi-lo)-.55)/.45,0,1);w=float(6*t**5-15*t**4+10*t**3);s['height_mm']=float(h+local[2]*w);s['boundary_local_mm']=(np.array(s['boundary_local_mm'])+local[:2]*w).tolist();s['pole_local_mm']=(np.array(s['pole_local_mm'])+local[:2]*w).tolist();changes.append({'section':i,'blend_weight':w,'displacement_display_mm':(delta*w).tolist()})
    newheights=[r['height_mm'] for r in profile['sections']];assert np.all(np.diff(newheights)>0);profile['body_height_range_mm']=[newheights[0],newheights[-1]];profile['body_axial_span_mm']=newheights[-1]-newheights[0];profile['engineering_revision']={'scope':'new neutral muscle-end geometry for direct forearm attachment, not elastic deformation result','source_profile_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_search_sha256':hashlib.sha256((OUT/'direct_muscle_interfaces/flexor_attachment_search.json').read_bytes()).hexdigest(),'distal_delta_display_mm':delta.tolist(),'proximal_boundary_unchanged':True,'bone_centres_changed':False,'sections':changes};(DEST/'FJ1486_profile.json').write_text(json.dumps(profile,indent=2),encoding='utf8');print('direct flexor distal change mm',delta,'norm',np.linalg.norm(delta),flush=True)

if __name__=='__main__':run()
