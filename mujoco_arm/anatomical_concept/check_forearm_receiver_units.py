"""Finite local gimbal articulation, shoulder seating and fitting checks."""
from build_forearm_receivers import *

def run():
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;d=json.loads((folder/'design.json').read_text());offset=d['standoff_mm'];parts={k:load(folder/(k+'.obj'))[0] for k in d['parts']};motions=[]
        for label,poses in [('yaw_revolution',[(a,0) for a in range(361)]),('pitch_at_5_yaws',[(yaw,pitch) for yaw in [0,45,90,180,270] for pitch in range(-80,81)])]:
            first=None;maxoverlap=0
            for yaw,pitch in poses:
                moved=articulate(parts,yaw,pitch,offset);h=[]
                for (a,sa),(b,sb) in itertools.combinations(moved.items(),2):
                    # Invariant same-stage pairs checked at zero. All pairs
                    # between different stages remain in every sampled pose.
                    stage=lambda k:2 if k=='pitch_crimp_blank' else 1 if k in YAW_PARTS else 0
                    if stage(a)==stage(b) and (yaw,pitch)!=(0,0):continue
                    v=intersection(sa,sb)
                    if v>.001:h.append({'pair':[a,b],'volume_mm3':v});maxoverlap=max(maxoverlap,v)
                if h and first is None:first={'yaw_deg':yaw,'pitch_deg':pitch,'hits':h}
            motions.append({'path':label,'samples':len(poses),'angular_step_deg':1,'return_same_rigid_configurations':True,'first_limit':first,'max_overlap_mm3':maxoverlap})
        h=5.5+offset
        checks=[]
        for label,body,seat in [('yaw_rear_shoulder','yaw_sleeve','yaw_base_washer'),('yaw_front_shoulder','yaw_sleeve','yaw_front_washer'),('pitch_shoulder_on_fork','pitch_sleeve','rotor_yoke')]:
            gap=float(parts[body].min_gap(parts[seat],1));checks.append({'label':label,'gap_mm_capped1':gap,'nominal_contact':gap<.0001})
        # Installation order on bench: loose pitch sleeve inside crimp eye,
        # both enter open yoke from +X; transverse pin follows, yaw components
        # assembled onto lever only after the yoke is ready.
        fixed=['rotor_yoke'];paths=[]
        for label,keys,direction,length in [('pitch_hub_and_sleeve_into_open_yoke',['pitch_crimp_blank','pitch_sleeve'],np.array([1,0,0]),18),('pitch_back_washer_and_loose_nut',['pitch_back_washer','pitch_nut'],np.array([0,-1,0]),12),('pitch_pin_and_front_washer',['pitch_M2x12','pitch_front_washer'],np.array([0,1,0]),15)]:
            first=None
            for distance in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
                hits=[]
                for key in keys:
                    for other in fixed:
                        v=intersection(parts[key].translate(direction*distance),parts[other])
                        if v>.001:hits.append({'moving':key,'fixed':other,'volume_mm3':v})
                if hits and first is None:first={'distance_mm':float(distance),'hits':hits}
            paths.append({'label':label,'moving':keys,'direction':direction.tolist(),'first_limit':first,'final_5mm_step_mm':.1,'other_step_mm':1});fixed+=keys
        report={'scope':'local rigid yaw/pitch sampling,shoulder contact and3 bench insertion paths only; no installed arm/route,thread,friction,preload or load proof','motion_paths':motions,'nominal_seat_contacts':checks,'assembly_paths':paths,'all_local_motion_clear':not any(r['first_limit'] for r in motions),'all_local_fitting_clear':not any(r['first_limit'] for r in paths),'source_sha256':{k:hashlib.sha256((folder/(k+'.obj')).read_bytes()).hexdigest() for k in parts},'frozen':False};(folder/'unit_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(fid,'local motion',report['all_local_motion_clear'],'fitting',report['all_local_fitting_clear'],flush=True)

if __name__=='__main__':run()
