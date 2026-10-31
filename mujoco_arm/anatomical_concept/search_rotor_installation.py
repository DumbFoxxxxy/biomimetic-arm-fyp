"""Search actual bare-elbow insertion paths without modifying either bone."""
from check_rotor_pocket import ROOT,DEST
from build_elbow_muscle_mounts import OUT,load,np,unit
from scipy.spatial.transform import Rotation
import json,hashlib

def run():
    base=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
    data=json.loads((base/'assembly.json').read_text(encoding='utf8'));axis=np.array(data['axis_source'])
    from build_reference import SHIFT
    E=np.array(data['centre_source_mm'])-SHIFT
    hum=load(DEST/'humerus_rotor_pocket.obj')[0];fore=load(DEST/'forearm_with_rotor.obj')[0];bush=load(base/'elbow_bearing_bush.obj')[0]
    prior=json.loads((base/'assembly_and_retention_checks.json').read_text(encoding='utf8'));side=np.array(prior['chosen_cartridge_path']['direction']);up=np.cross(axis,side)
    rows=[];chosen=None
    def contact(parts,delta):
        hits=[]
        for name,part in parts.items():
            volume=max(0,float((part.translate(delta)^hum).volume()))
            if volume>.001:hits.append({'part':name,'volume_mm3':volume})
        return hits
    # Prior tested direction first; try other orientations in the elbow-normal
    # plane, then skewed axial approaches. No adjacent collision pair is masked.
    for skew in [0,.5,-.5,1.,-1.]:
        if chosen:break
        for angle in [0,90,150,45,120]:
            if chosen:break
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix()
            parts={'forearm_with_rotor':fore.transform(np.c_[R,E-R@E]),'bush':bush.transform(np.c_[R,E-R@E])}
            for azimuth in [0,180,90,270,30,330,60,300,120,240,150,210,15,345,45,315,75,285,105,255,135,225,165,195]:
                direction=unit(side*np.cos(np.radians(azimuth))+up*np.sin(np.radians(azimuth))+axis*skew)
                samples=[];clear=True
                for offset in range(60,-1,-1):
                    hits=contact(parts,direction*offset);samples.append({'offset_mm':offset,'hits':hits})
                    if hits:clear=False;break
                row={'elbow_flexion_deg':angle,'azimuth_deg':azimuth,'axial_skew':skew,'direction':direction.tolist(),'coarse_step_mm':1,'first_interference':next((r for r in samples if r['hits']),None),'coarse_clear':clear}
                rows.append(row)
                if clear:
                    fine=[{'offset_mm':float(offset),'hits':contact(parts,direction*offset)} for offset in np.linspace(5,0,51)]
                    if not any(r['hits'] for r in fine):chosen=dict(row,samples=samples,fine_samples=fine,last_5mm_step_mm=.1);break
            print('search',len(rows),'skew',skew,'flexion',angle,'found',chosen is not None,flush=True)
    result={'scope':'bare elbow forearm+rotor/bush straight translation search with unchanged humerus; sleeves/pin absent. Local insertion only, not complete assembly or driver proof.','trials':rows,'chosen':chosen,'geometry_sha256':{n:hashlib.sha256((DEST/n).read_bytes()).hexdigest() for n in ['humerus_rotor_pocket.obj','forearm_with_rotor.obj']},'frozen':False}
    (DEST/'alternate_installation_search.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print('chosen',None if chosen is None else {k:v for k,v in chosen.items() if k not in ['samples','fine_samples']},flush=True)
if __name__=='__main__':run()
