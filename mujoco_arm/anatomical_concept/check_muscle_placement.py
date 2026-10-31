"""Neutral context only: hard collisions block placement, TPU contacts stay unknown."""
from parameterize_elbow_muscles import *
from build_AC_connection import load
SKELETON=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'

def run():
    manifest=json.loads((SKELETON/'whole_arm_STEP_check.json').read_text(encoding='utf8'));bones={};boxes={}
    for key,record in manifest['parts'].items():
        path=ROOT/record['source_STEP'];obj=path.parent.parent/(record['part']+'.obj');s=load(obj)[0];bones[key]=s;boxes[key]=np.array(s.bounding_box())
    results=[]
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;data=json.loads((folder/'design.json').read_text(encoding='utf8'));parts={name:load(folder/(name+'.obj'))[0] for name in data['components']}
        def overlaps(offset,soft=False):
            hits=[]
            for name,s in parts.items():
                if (name=='TPU_muscle')!=soft:continue
                moved=s.translate(offset);box=np.array(moved.bounding_box())
                for key,body in bones.items():
                    delta=np.maximum(box[:3]-boxes[key][3:],boxes[key][:3]-box[3:])
                    if np.any(delta>=0):continue
                    volume=max(0,float((moved^body).volume()))
                    if volume>.001:hits.append({'module_part':name,'skeleton_part':key,'volume_mm3':volume})
            return hits
        neutral=overlaps(np.zeros(3));soft=overlaps(np.zeros(3),True);direction=np.array([0.,-1. if fid=='FJ1486' else 1.,0.]);offsets=[];first=None
        for distance in range(0,33,2):
            hits=overlaps(direction*distance);offsets.append({'distance_mm':distance,'hard_interferences':hits})
            if not hits:first=distance;break
        result={'id':fid,'scope':'unmounted common-coordinate neutral placement screen; no hardware attachment, motion or TPU deformation pass','neutral_hard_interferences':neutral,'neutral_TPU_intersections':soft,'soft_interpretation':'intersection volume is a contact/deformation review flag, neither automatic failure nor acceptance','outward_translation_direction':direction.tolist(),'tested_offsets':offsets,'first_sampled_rigid_clear_offset_mm':first,'offset_adopted':False,'why_not_adopted':'bone adapters, tendon route, soft contacts and full motion still unvalidated','frozen':False}
        results.append(result);print(fid,'neutral hard overlaps',len(neutral),'soft flags',len(soft),'first rigid-clear offset',first,flush=True)
    (DEST/'neutral_placement_screen.json').write_text(json.dumps({'scope':'neutral placement diagnostics only','muscles':results},indent=2),encoding='utf8')

if __name__=='__main__':run()
