"""Undeformed BCC-inspired PSM coupon; no claim of verified TPU contraction."""
from pathlib import Path
import sys,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tmp/shoulder_motion_design_20261029/vendor'),str(ROOT/'mujoco_arm')]
import manifold3d as md
from stl_io import write_binary_stl
OUT=ROOT/'mujoco_arm/exports/tpu_muscle_design'

def cylinder(a,b,r):
    a=np.array(a,dtype=float);b=np.array(b,dtype=float);axis=b-a;length=np.linalg.norm(axis);axis/=length
    helper=np.array([0.,0.,1.]) if abs(axis[2])<.9 else np.array([0.,1.,0.]);u=np.cross(helper,axis);u/=np.linalg.norm(u);v=np.cross(axis,u)
    return md.Manifold.cylinder(length,r,circular_segments=24,center=True).transform(np.column_stack([np.column_stack([u,v,axis]),(a+b)/2]))

def main():
    OUT.mkdir(parents=True,exist_ok=True);x,y,z,p=15.,10.,30.,2.;count=5;parts=[]
    for i in range(count):
        centre=np.array([(i+.5)*x,0.,0.])
        for xx in [i*x,(i+1)*x]:
            for yy in [-y/2,y/2]:
                for zz in [-z/2,z/2]:
                    corner=np.array([xx,yy,zz]);parts.append(cylinder(centre,corner,p/2));parts.append(md.Manifold.sphere(p/2,circular_segments=16).translate(corner))
        parts.append(cylinder(centre-[1.25,0,0],centre+[1.25,0,0],2.2))
    for start in [-10.,count*x]:parts.append(md.Manifold.cube((10.,y+2,z+2)).translate([start,-y/2-1,-z/2-1]))
    solid=md.Manifold.batch_boolean(parts,md.OpType.Add)-cylinder([-11,0,0],[count*x+11,0,0],.9)
    for start in [-10.,count*x]:
        for zz in [-10.,10.]:solid=solid-cylinder([start-1,0,zz],[start+11,0,zz],1.6)
    mesh=solid.to_mesh();vertices=np.asarray(mesh.vert_properties)[:,:3];faces=np.asarray(mesh.tri_verts)
    # This STL intentionally uses millimetres, unlike the old MuJoCo assets.
    write_binary_stl(OUT/'psm_reference_coupon_mm.stl',vertices[faces],'PSM-inspired undeformed coupon; MILLIMETRES; not contraction-validated')
    report={'scope':'Undeformed printable geometry candidate, not an exact reproduction of author CAD and not a verified muscle actuator',
        'topology':'Figure-5-inspired body-centre to eight corners, repeated along X; discrete bored wire guide rings; dense end blocks',
        'units':'mm','cell_xyz_mm':[x,y,z],'branch_diameter_mm':p,'cell_count_x':count,'active_pitch_length_mm':count*x,'end_blocks_each_mm':10.,'wire_guide_bore_diameter_mm':1.8,'mount_hole_diameter_mm':3.2,
        'material_proposed':'TPU; paper uses 95A, actual filament and printing process not selected','manifold_status':str(solid.status()),'positive_material_components':sum(s.volume()>1e-6 for s in solid.decompose()),'axial_bore_ray_hits':len(solid.ray_cast([-11,0,0],[count*x+11,0,0])),
        'volume_mm3':solid.volume(),'bounds_mm':list(solid.bounding_box()),'contraction_verified':False,'force_model_calibrated':False,
        'paper_psm_2020_example':{'overall_length_mm':250,'contracted_length_mm':115,'reported_overall_shortening_fraction':.54,'not_a_rating_for_this_coupon':True},
        'twenty_percent_note':'AIS2026 calls ~20% a biological-muscle macroscopic reference; it is not a universal PSM contraction ceiling.'}
    if report['manifold_status']!='Error.NoError' or report['positive_material_components']!=1 or report['axial_bore_ray_hits']!=0:raise RuntimeError(report)
    (OUT/'psm_reference_coupon.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    source=ROOT/'mujoco_arm/exports/blade_motion_design/full_arm_cable_coordinates.csv';rows=[]
    for row in csv.DictReader(source.open(encoding='utf-8-sig')):
        stroke=float(row['required_geometric_excursion_mm']);out={'reference_cable_id':row['id'],'reference_wire_stroke_mm':stroke}
        for eps in [.18,.20,.35,.50]:out[f'min_active_length_at_{int(eps*100)}pct_mm']=stroke/eps
        if row['id'].startswith('S'):
            # Paper Table I reports both shortening AND extension about its
            # nominal whole-muscle length. Do not force Lmax to mean Lrest.
            out['paper_scaled_nominal_whole_length_lower_mm']=float(row['maximum_active_length_mm'])/(380./250.)
            out['paper_scaled_nominal_whole_length_upper_mm']=float(row['minimum_active_length_mm'])/(115./250.)
        else:
            out['paper_scaled_nominal_whole_length_lower_mm']='not applicable: wire includes drum arc'
            out['paper_scaled_nominal_whole_length_upper_mm']='not applicable: wire includes drum arc'
        out['status']='Screening only: assumes 1:1 muscle shortening vs active wire stroke; old anchors not approved for TPU'
        rows.append(out)
    with (OUT/'reference_stroke_budget.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('STROKE_SCREEN',[(r['reference_cable_id'],round(r['reference_wire_stroke_mm'],2),round(r['min_active_length_at_20pct_mm'],1)) for r in rows],flush=True)
if __name__=='__main__':main()
