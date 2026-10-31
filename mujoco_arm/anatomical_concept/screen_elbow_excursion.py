"""Direct-anchor diagnostic and ideal-drum stroke requirements, not a drive pass."""
from build_elbow_connection import load,unit,rod
from parameterize_elbow_muscles import ROOT,OUT,DEST,SHIFT
from extract_attachment_regions import projector
from scipy.spatial.transform import Rotation
import numpy as np,json
SKELETON=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'

def run():
    definition=json.loads((OUT/'joint_model/definition.json').read_text(encoding='utf8'));a=json.loads((SKELETON/'assembly.json').read_text(encoding='utf8'))
    E=np.array(a['centre_source_mm'])-SHIFT;axis=np.array(a['axis_source']);up=unit(np.array(definition['centres_source_mm']['GH'])-SHIFT-E)
    parts={name:load(SKELETON/(name+'.obj')) for name in a['components']};rows=[]
    for fid in ['FJ1486','FJ1477']:
        anchors=[]
        for role,structure,target in [('origin','humerus','humerus_GH_elbow_fork'),('insertion','ulna','forearm_elbow_hub')]:
            patch=np.load(OUT/'attachment_regions'/(fid+'_'+role+'_'+structure+'.npz'))['bone_points_source_mm']-SHIFT
            if role=='origin':
                height=(patch-E)@up;patch=patch[height>=np.percentile(height,90)]
            mean=patch.mean(0);s,v,f=parts[target];point,d,indices=projector(v,f)(mean[None]);i=int(indices[0]);tri=v[f[i]];normal=unit(np.cross(tri[1]-tri[0],tri[2]-tri[0]));anchors.append(point[0]+normal*3)
        origin,insertion=anchors;lengths=[];interferences=[]
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();end=E+R@(insertion-E);length=np.linalg.norm(end-origin);direction=(end-origin)/length
            # Endpoint 5mm neighborhoods are diagnostic terminal regions, not
            # whole neighboring bones excluded from collision checks.
            steel=rod(origin+direction*5,end-direction*5,.3);hits=[]
            for name,(s,v,f) in parts.items():
                if a['part_owners'][name]=='elbow':s=s.transform(np.c_[R,E-R@E])
                volume=max(0,float((steel^s).volume()))
                if volume>.001:hits.append({'part':name,'volume_mm3':volume})
            lengths.append(float(length));interferences.append({'angle_deg':angle,'hits':hits})
        moment=-np.gradient(lengths,np.radians(np.arange(151)));desired=1 if fid=='FJ1486' else -1;wrong=np.where(moment*desired<=0)[0]
        profile=json.loads((DEST/(fid+'_profile.json')).read_text(encoding='utf8'));span=profile['body_axial_span_mm']
        row={'id':fid,'anchor_scope':'illustrative atlas/proximity-based engineering probes, NOT installed or frozen interfaces','origin_probe_display_mm':origin.tolist(),'insertion_probe_display_mm':insertion.tolist(),'standoff_mm':3,'direct_route_length_mm':lengths,'direct_route_geometric_excursion_mm':max(lengths)-min(lengths),'moment_arm_negative_dL_dtheta_mm':moment.tolist(),'first_wrong_or_zero_torque_sign_deg':int(wrong[0]) if len(wrong) else None,'colliding_samples':sum(bool(p['hits']) for p in interferences),'first_colliding_sample':next((p for p in interferences if p['hits']),None),'interferences':interferences,'source_profile_body_span_mm':span,'ideal_drum_comparisons':[{'radius_mm':r,'rotation_deg':150,'required_linear_stroke_mm':float(r*np.radians(150)),'stroke_fraction_of_profile_span':float(r*np.radians(150)/span),'physical_drum_exists':False,'TPU_capability_verified':False} for r in [8,10,12,14]],'drive_pass':False}
        rows.append(row);print(fid,'direct excursion',round(row['direct_route_geometric_excursion_mm'],2),'wrong sign',row['first_wrong_or_zero_torque_sign_deg'],'colliding poses',row['colliding_samples'],flush=True)
    report={'scope':'geometry diagnostic only; straight line is not a validated shaped muscle/fiber path. Ideal drums are comparisons, no nonexistent guide or drum imposes a constraint.','angular_step_deg':1,'moment_arm_method':'finite differences of direct endpoint distance in radians, not force or compliance simulation','wire_diameter_mm':.6,'endpoint_neighborhood_trim_mm':5,'collision_volume_threshold_mm3':.001,'muscles':rows,'frozen':False}
    (DEST/'excursion_screen.json').write_text(json.dumps(report,indent=2),encoding='utf8')

if __name__=='__main__':run()
