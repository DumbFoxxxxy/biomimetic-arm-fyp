"""Direct passive tendon geometry and assumed axial muscle-end travel.

Both endpoint holes are real M4 geometry. Axial cap translation is a diagnostic
assumption, NOT an imposed guide, demonstrated TPU deformation or drive pass.
"""
from build_elbow_side_lever import *
from scipy.spatial.transform import Rotation

def run():
    mounts=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'));offsets={p['id']:np.array(p['translation_display_mm']) for p in mounts['placements']};modules={}
    for fid in ['FJ1486','FJ1477']:
        data=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));F=np.array(data['basis_columns_display']);pad=np.array(data['pad_centres_display_mm'][1])+offsets[fid]
        modules[fid]={'frame':F,'eyelet_centre':pad+F@np.array([0,4,4.15])}
    results=[]
    for variant in json.loads((DEST/'variants.json').read_text(encoding='utf8'))['variants']:
        folder=DEST/variant;lever=json.loads((folder/'design.json').read_text(encoding='utf8'));E=np.array(lever['joint_centre_display_mm']);axis=np.array(lever['joint_axis']);rows=[]
        for fid,module in modules.items():
            F=module['frame'];z=F[:,2]
            for eyelet,anchor in lever['eyelets_display_mm'].items():
                B0=np.array(anchor)-axis*np.sign(lever['eyelet_axis_station_mm'])*lever.get('free_span_entry_offset_towards_cap_mm',2.4);side=1 if (B0-module['eyelet_centre'])@F[:,0]>=0 else -1
                # Bore is along cap-localX. Exit lies outside11mm cap radius;
                # constant interior tendon span is independent of axial travel.
                P0=module['eyelet_centre']+F[:,0]*side*11.3;length=float(np.linalg.norm(B0-P0));path=[];last=0.;invalid=None
                for angle in range(151):
                    R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();B=E+R@(B0-E);D=B-P0;dot=float(D@z);discriminant=dot*dot-float(D@D)+length*length
                    candidates=[] if discriminant<0 else [-dot+np.sqrt(max(0,discriminant)),-dot-np.sqrt(max(0,discriminant))]
                    if not candidates:
                        invalid=angle;break
                    stroke=min(candidates,key=lambda s:abs(s-last));last=float(stroke);P=P0-z*stroke
                    path.append({'angle_deg':angle,'cap_axial_shortening_mm':float(stroke),'cap_exit_display_mm':P.tolist(),'lever_anchor_display_mm':B.tolist(),'passive_free_span_mm':float(np.linalg.norm(B-P)),'quadratic_discriminant_mm2':discriminant})
                strokes=np.array([p['cap_axial_shortening_mm'] for p in path]);derivative=np.gradient(strokes,np.radians(np.arange(len(path)))) if len(path)>1 else np.array([0.])
                desired=1 if fid=='FJ1486' else -1;wrong=np.flatnonzero(derivative*desired<=0)
                # This derivative is virtual axial stroke per radian only; no
                # anatomical fiber or real muscle mechanical advantage claim.
                row={'muscle':fid,'lever_eyelet':eyelet,'cap_eyelet_exit_side':side,'neutral_free_tendon_span_mm':length,'path':path,'first_no_axial_solution_deg':invalid,'first_wrong_or_zero_desired_stroke_sign_deg':int(wrong[0]) if len(wrong) else None,'axial_shortening_range_mm':[float(strokes.min()),float(strokes.max())] if len(strokes) else None,'axial_excursion_mm':float(np.ptp(strokes)) if len(strokes) else None,'axial_stroke_per_radian_mm':derivative.tolist(),'assumed_path_clearance_checked':False,'TPU_deformation_verified':False,'complete_sign_consistent_axial_family':len(path)==151 and not len(wrong),'drive_pass':False}
                rows.append(row)
                print(variant,fid,eyelet,'samples',len(path),'wrong',row['first_wrong_or_zero_desired_stroke_sign_deg'],'stroke',row['axial_shortening_range_mm'],flush=True)
        report={'variant':variant,'scope':'physical hole references with an UNVALIDATED axial cap path family; no nonexistent guide is imposed','equation':'|B(theta)-(P0-z*s)|=|B(0)-P0|; nearest continuous quadratic root from s=0','modules':rows,'source_lever_design_sha256':hashlib.sha256((folder/'design.json').read_bytes()).hexdigest(),'frozen':False}
        (folder/'tendon_axial_diagnostic.json').write_text(json.dumps(report,indent=2),encoding='utf8');results.append({'variant':variant,'sign_consistent_combinations':[{'muscle':r['muscle'],'eyelet':r['lever_eyelet'],'excursion_mm':r['axial_excursion_mm']} for r in rows if r['complete_sign_consistent_axial_family']]})
    (DEST/'tendon_comparison.json').write_text(json.dumps({'scope':'assumed endpoint motion comparison, not flexible muscle or cable validation','variants':results},indent=2),encoding='utf8')
if __name__=='__main__':run()
