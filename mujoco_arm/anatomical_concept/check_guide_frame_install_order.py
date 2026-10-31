"""Alternative real installation order after front-bolt approach was blocked.

Install guides before forearm, then explicitly test forearm, spacers, pin and
lever paths with actual frame solids still present. Muscle installation remains
separate and unverified; nothing is deleted from the final geometry.
"""
from build_guide_ground_frames import *

def run():
    a=json.loads((DEST7/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};selection=json.loads((DEST/'paired_placement/selection.json').read_text())['placements'];results=[]
    guide_names={fid+'_'+k for fid in selection for k in json.loads((TEMPLATE/'design.json').read_text())['parts']};mount_names={fid+'_'+k for fid in selection for k in local_receiver()[1]}
    adapted={fid+'_bone_adapter' for fid in selection};bone_hardware={k for k in rows if k.startswith(('FJ1486_bone_','FJ1477_bone_'))}
    structural=[k for k in world if not k.startswith(('FJ1486_','FJ1477_')) or k in guide_names|mount_names|adapted|bone_hardware]
    elbow_pin=['elbow_spacer_left','elbow_spacer_right','elbow_shoulder_pin','elbow_washer_left','elbow_washer_right','elbow_retaining_nut'];elbow_moving=[k for k,r in rows.items() if r['owner']=='elbow']+['flexor_tendon_stop','extensor_tendon_stop']
    fixed=[k for k in structural if k not in elbow_moving+elbow_pin]
    def probe(label,moving,against,direction,length):
        first=None;count=0
        for offset in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
            h=[];count+=1
            for k in moving:h.extend({'moving':k,**hit} for hit in hits(world[k].translate(np.array(direction)*offset),{j:world[j] for j in against}))
            if h:first={'offset_mm':float(offset),'hits':h};break
        result={'label':label,'moving':moving,'fixed':against,'direction':np.asarray(direction).tolist(),'distance_mm':length,'sample_count':count,'first_obstruction':first,'clear':first is None,'final_5mm_step_mm':.1,'other_step_mm':1};results.append(result);print(label,'clear',first is None,first or '',flush=True)
    # Changed guide fastener paths before forearm exists; full mounted rigid
    # guide context retained. This is a stated sequence, not a collision exclusion.
    for fid,s in selection.items():
        F=np.array(s['frame_columns']);pair=[fid+'_key_M3x20',fid+'_key_front_washer'];probe(fid+'_key_front_bolt_BEFORE_forearm',pair,[k for k in fixed if k not in pair],F[:,2],25)
    skeleton=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection';old=json.loads((skeleton/'assembly_and_retention_checks.json').read_text());axis=np.array(json.loads((skeleton/'assembly.json').read_text())['axis_source'])
    pair=['FJ3391','elbow_bearing_bush'];probe('forearm_and_bush_AFTER_grounded_guides',pair,fixed,np.array(old['chosen_cartridge_path']['direction']),45);fixed+=pair
    for p in old['separate_sleeve_insertions']:
        key=p['part'];probe(key,[key],fixed,-np.array(p['direction']),85);fixed.append(key)
    pair=['elbow_shoulder_pin','elbow_washer_left'];probe('elbow_pin_after_frame',pair,fixed,-axis,old['pin_insertion']['distance_mm']);fixed+=pair
    for key in ['elbow_washer_right','elbow_retaining_nut']:probe(key,[key],fixed,axis,25);fixed.append(key)
    lever=json.loads((BASE/'design.json').read_text());n=np.array(lever['bone_anchors'][0]['bolt_axis']);probe('removable_lever_AFTER_grounded_guides',['side_lever'],fixed,n,60);fixed.append('side_lever')
    for i in [0,1]:
        prefix='bone_'+str(i);pair=[prefix+'_bolt',prefix+'_front_washer'];probe(prefix+'_front',pair,fixed,n,45);fixed+=pair
        pair=[prefix+'_nut',prefix+'_back_washer'];probe(prefix+'_rear',pair,fixed,-n,35);fixed+=pair
    report={'scope':'alternative sampled rigid sequence:grounded guide/key hardware BEFORE forearm, then actual forearm/cartridge, spacers,pin and lever; final muscles/cables are not removed from final assembly, but their subsequent installation and compliant handling are unverified','paths':results,'all_clear':all(r['clear'] for r in results),'complete_arm_assembly_verified':False,'source_sha256':{str(ROOT/r['file']):digest(ROOT/r['file']) for r in rows.values()},'frozen':False}
    (DEST7/'alternative_install_order.json').write_text(json.dumps(report,indent=2),encoding='utf8')

if __name__=='__main__':run()
