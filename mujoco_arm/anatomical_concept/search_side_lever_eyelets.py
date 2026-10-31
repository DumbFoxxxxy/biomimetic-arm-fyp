"""Search angular eyelet positions under an explicitly assumed axial cap path.
Outputs are candidate positions, not verified muscle deformation or guide laws.
"""
from build_elbow_side_lever import *
from scipy.spatial.transform import Rotation

def run():
    mounts=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'));offsets={p['id']:np.array(p['translation_display_mm']) for p in mounts['placements']};modules={}
    for fid in ['FJ1486','FJ1477']:
        data=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));F=np.array(data['basis_columns_display']);pad=np.array(data['pad_centres_display_mm'][1])+offsets[fid]
        modules[fid]=(F,pad+F@np.array([0,4,4.15]))
    rows=[]
    for variant in json.loads((DEST/'variants.json').read_text(encoding='utf8'))['variants']:
        folder=DEST/variant;lever=json.loads((folder/'design.json').read_text(encoding='utf8'));E=np.array(lever['joint_centre_display_mm']);axis=np.array(lever['joint_axis']);down=np.array(lever['zero_forearm_direction']);radius=lever['nominal_tendon_reference_radius_mm'];station=lever['eyelet_axis_station_mm'];candidates=[]
        for fid,(F,centre) in modules.items():
            found=[];desired=1 if fid=='FJ1486' else -1
            for phase in range(360):
                R=Rotation.from_rotvec(np.radians(np.arange(151)+phase)[:,None]*axis).as_matrix();hole=E+axis*station+np.einsum('nij,j->ni',R,down*radius);B=hole-axis*np.sign(station)*2.4
                side=1 if (B[0]-centre)@F[:,0]>=0 else -1;P0=centre+F[:,0]*side*11.3;D=B-P0;dot=D@F[:,2];L=np.linalg.norm(D[0]);disc=dot*dot-np.sum(D*D,axis=1)+L*L
                if disc.min()<0:continue
                branch=1 if dot[0]>=0 else -1;s=-dot+branch*np.sqrt(np.maximum(disc,0));derivative=np.gradient(s,np.radians(np.arange(151)))*desired
                if derivative.min()<=1e-5:continue
                found.append({'phase_deg':phase,'eyelet_display_mm':hole[0].tolist(),'free_span_entry_display_mm':B[0].tolist(),'axial_excursion_mm':float(np.ptp(s)),'minimum_desired_axial_stroke_per_radian_mm':float(derivative.min()),'minimum_quadratic_discriminant_mm2':float(disc.min()),'neutral_free_tendon_span_mm':float(L),'cap_exit_side':side,'axial_shortening_range_mm':[float(s.min()),float(s.max())]})
            found.sort(key=lambda p:(-p['minimum_desired_axial_stroke_per_radian_mm'],p['axial_excursion_mm']))
            chosen=found[0] if found else None
            candidates.append({'muscle':fid,'feasible_phase_samples':len(found),'selected_diagnostic_candidate':chosen,'all_sign_consistent_candidates':found})
            print(variant,fid,'feasible phases',len(found),'best',chosen,flush=True)
        rows.append({'variant':variant,'muscles':candidates})
    (DEST/'eyelet_phase_search.json').write_text(json.dumps({'scope':'diagnostic angular search under unvalidated axial cap translation; no guide is imposed','angular_phase_step_deg':1,'elbow_path_step_deg':1,'selection':'maximize minimum desired axial stroke derivative over0..150deg; then smaller excursion. No force/stress/muscle feasibility claim.','variants':rows,'frozen':False},indent=2),encoding='utf8')
if __name__=='__main__':run()
