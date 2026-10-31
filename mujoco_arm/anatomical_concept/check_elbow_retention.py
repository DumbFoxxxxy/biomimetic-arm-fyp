"""Finite assembly-path and retention probes of exported E1 parts."""
from build_elbow_connection import *

def run():
    data=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));parts={n:load(DEST/(n+'.obj'))[0] for n in data['components']}
    axis=np.array(data['axis_source']);E=np.array(data['centre_source_mm'])-SHIFT
    definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'));up=unit(np.array(definition['centres_source_mm']['GH'])-SHIFT-E);up=unit(up-axis*(up@axis));side=np.cross(axis,up)
    def hits(moving,fixed,delta):
        result=[]
        for a in moving:
            moved=parts[a].translate(delta)
            for b in fixed:
                volume=max(0,float((moved^parts[b]).volume()))
                if volume>.001:result.append({'pair':[a,b],'volume_mm3':volume})
        return result
    forearm=['forearm_elbow_hub','elbow_bearing_bush'];fixed=[n for n in parts if n not in forearm];rows=[]
    for label,direction in [('axial_plus',axis),('axial_minus',-axis),('up',up),('down',-up),('side',side),('minus_side',-side)]:
        first=None
        for distance in np.arange(.1,2.001,.1):
            contact=hits(forearm,fixed,direction*distance)
            if contact:first={'translation_mm':float(distance),'hits':contact};break
        rows.append({'direction':label,'first_interference':first})
    pin_bundle=['elbow_shoulder_pin','elbow_retaining_nut','elbow_washer_left','elbow_washer_right'];pin_rows=[]
    for sign in [-1,1]:
        first=None
        for distance in np.arange(.1,2.001,.1):
            contact=hits(pin_bundle,['humerus_GH_elbow_fork'],axis*sign*distance)
            if contact:first={'translation_mm':float(distance),'hits':contact};break
        pin_rows.append({'axis_sign':sign,'first_interference':first})
    # Insert shank with its head-side washer before fitting the nut/end washer.
    insertion=[];travel=int(np.ceil(data['custom_pin_shank_length_mm']+10))
    fixed_insert=[n for n in parts if n not in pin_bundle]
    for distance in range(travel,-1,-1):
        contact=hits(['elbow_shoulder_pin','elbow_washer_left'],fixed_insert,-axis*distance)
        insertion.append({'offset_mm':distance,'hits':contact})
    # Radial installation of forearm/bush/spacers with the pin absent.
    direction=unit(np.array(data['forearm_bridges_display_mm'][0]).mean(0)-E)
    cartridge=forearm+['elbow_spacer_left','elbow_spacer_right'];cartridge_path=[]
    for distance in range(45,-1,-1):
        cartridge_path.append({'offset_mm':distance,'hits':hits(cartridge,['humerus_GH_elbow_fork'],direction*distance)})
    search=[];chosen=None
    for flexion in [0,90,120,150,45,60]:
        if chosen:break
        R=Rotation.from_rotvec(axis*np.radians(flexion)).as_matrix()
        transformed={n:parts[n].transform(np.c_[R,E-R@E]) for n in forearm}
        for azimuth in [0,90,180,270,45,135,225,315,15,30,60,75,105,120,150,165,195,210,240,255,285,300,330,345]:
            direction=up*np.cos(np.radians(azimuth))+side*np.sin(np.radians(azimuth));first=None;samples=[]
            for distance in range(45,-1,-1):
                contact=[]
                for name,s in transformed.items():
                    volume=max(0,float((s.translate(direction*distance)^parts['humerus_GH_elbow_fork']).volume()))
                    if volume>.001:contact.append({'part':name,'volume_mm3':volume})
                samples.append({'offset_mm':distance,'hits':contact})
                if contact:first=samples[-1];break
            search.append({'elbow_flexion_deg':flexion,'insertion_direction':direction.tolist(),'first_interference':first})
            if first is None:
                # Refine the last 5mm; keep the coarse approach samples as evidence.
                fine=[]
                for distance in np.linspace(5,0,51):
                    contact=[]
                    for name,s in transformed.items():
                        volume=max(0,float((s.translate(direction*distance)^parts['humerus_GH_elbow_fork']).volume()))
                        if volume>.001:contact.append({'part':name,'volume_mm3':volume})
                    fine.append({'offset_mm':float(distance),'hits':contact})
                if any(s['hits'] for s in fine):continue
                chosen={'moving_parts':forearm,'elbow_flexion_deg':flexion,'direction':direction.tolist(),'step_mm':1,'final_5mm_step_mm':.1,'samples':samples,'fine_samples':fine};break
    sleeves=[]
    for name,sign in [('elbow_spacer_left',-1),('elbow_spacer_right',1)]:
        samples=[]
        fixed=[n for n in parts if n not in pin_bundle and n!=name]
        for distance in range(85,-1,-1):samples.append({'offset_mm':distance,'hits':hits([name],fixed,axis*sign*distance)})
        sleeves.append({'part':name,'direction':(-axis*sign).tolist(),'step_mm':1,'samples':samples,'collision_free':not any(s['hits'] for s in samples)})
    result={'scope':'finite straight geometric retention and ideal assembly paths; threads/preload, general escape, strength and tolerance not verified','translation_step_mm':.1,'forearm_retention':rows,'all_forearm_directions_blocked_within_2mm':all(r['first_interference'] is not None for r in rows),'pin_bundle_retention':pin_rows,'both_pin_directions_blocked_within_2mm':all(r['first_interference'] is not None for r in pin_rows),'pin_insertion':{'step_mm':1,'distance_mm':travel,'collision_free':not any(r['hits'] for r in insertion),'samples':insertion},'cartridge_insertion':{'step_mm':1,'distance_mm':45,'collision_free':not any(r['hits'] for r in cartridge_path),'samples':cartridge_path},'frozen':False}
    result.update(cartridge_direction_search=search,chosen_cartridge_path=chosen)
    result.update(separate_sleeve_insertions=sleeves,sequential_installation_geometry_pass=chosen is not None and all(s['collision_free'] for s in sleeves) and result['pin_insertion']['collision_free'],assembly_sequence=['forearm+bush via chosen radial path, pin and sleeves absent','slide two sleeves in axially through fork holes','insert shoulder pin with head-side washer','fit far washer and retaining nut; threads/preload not validated'])
    (DEST/'assembly_and_retention_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print('forearm retention',result['all_forearm_directions_blocked_within_2mm'],'pin capture',result['both_pin_directions_blocked_within_2mm'],'pin insertion',result['pin_insertion']['collision_free'],'cartridge insertion',result['cartridge_insertion']['collision_free'])
    print('chosen forearm insertion',None if chosen is None else {k:v for k,v in chosen.items() if k not in ['samples','fine_samples']},'tested',len(search),'sequential assembly',result['sequential_installation_geometry_pass'])

if __name__=='__main__':run()
