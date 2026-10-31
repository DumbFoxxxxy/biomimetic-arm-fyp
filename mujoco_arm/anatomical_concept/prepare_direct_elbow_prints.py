"""Independent bench-print candidates for the current direct muscle geometry.

Standalone distal caps are for bench specimens only. The installed arm uses
the single integrated yoke. No nozzle, support or slicing qualification here.
"""
from build_direct_muscle_yoke import *
import prepare_muscle_prints as bench

def run():
    assert ('--access-relief' in sys.argv or EXTENSOR_REANCHOR) and SHIFTED
    for fid in ['FJ1486','FJ1477']:
        bench.DEST=muscle_folder(fid).parent;bench.run((fid,))
    folder=DEST/'print_candidates';folder.mkdir(exist_ok=True)
    s,v,f=load(DEST/'direct_muscle_yoke.obj');_,axes=np.linalg.eigh(np.cov(v.T));R=axes.T
    if np.linalg.det(R)<0:R[0]*=-1
    transformed=v@R.T;t=np.r_[-(transformed.max(0)[:2]+transformed.min(0)[:2])/2,-transformed[:,2].min()]
    placed=s.transform(np.c_[R,t]);meta=write(folder/'direct_muscle_yoke.stl',placed)
    result={'scope':'print placement only, no slicer support or strength qualification','part':'direct_muscle_yoke','assembly_source':str((DEST/'direct_muscle_yoke.obj').relative_to(ROOT)),'assembly_source_sha256':hashlib.sha256((DEST/'direct_muscle_yoke.obj').read_bytes()).hexdigest(),'original_to_print_affine_mm':np.c_[R,t].tolist(),'bounds_mm':list(placed.bounding_box()),'print_mesh':meta,'requires_support_review':True,'material_selected':False,'installed_distal_caps':'integral with this yoke; do not install separate bench distal caps','frozen':False}
    (folder/'print_transforms.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    template={'scope':'blank physical test record; values must be measurements, not predicted travel','specimen_id':None,'muscle_id':None,'printer':None,'nozzle_mm':None,'TPU_brand_batch':None,'shore_A':None,'slicer_profile':None,'drying':None,'measured_strut_diameter_mm':None,'steel_diameter_mm':None,'conduit_reaction_clamped':False,'active_crimp_grip_tested':False,'preload_N':None,'trials':[],'trial_fields':['cycle','wire_feed_displacement_mm','end_to_end_distance_mm','elapsed_s','measured_wire_tension_N_or_null','return_residual_mm','wire_jam','end_damage','notes'],'all_measured':False}
    target=DEST/'physical_test_record_TEMPLATE.json'
    if not target.exists():target.write_text(json.dumps(template,indent=2),encoding='utf8')
    print('bench sets prepared; separate distal caps are not installed-arm parts',flush=True)

if __name__=='__main__':run()
