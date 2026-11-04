"""Map rotated dual-hole patterns against current anatomical bone geometry.

Only hole axes rotate: original cap centre, normal, steel port and bones stay
fixed. Rays establish geometric candidates, not complete mounts or strength.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,shoulder_pads

DEST=A.DESIGN/'s1_deltoid_origin_holes'


def main():
    DEST.mkdir(exist_ok=True)
    file=A.DESIGN/'s1_teres_assembly_order/assembly_check.json';parent=read(file)
    check=file.parent/'installation_checks.json';proof=read(check)
    assert proof['assembly_sha256']==sha(file) and proof['nominal_sequence_and_capture_geometry_pass']
    for row in parent['sources'].values():assert sha(ROOT/row['file'])==row['sha256']
    bones={key:load_component(parent['sources'][key]) for key in ['bone:scapula','bone:clavicle','bone:humerus']}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=read(route_file)['routes'];groups=[]
    for head,target in [('FJ1468','bone:clavicle'),('FJ1513','bone:scapula'),('FJ1467','bone:scapula')]:
        route=next(row for row in routes if row['head']==head);port=route['ports']['origin']
        I=np.array(port['inner_face_source_mm']);n=A.D.unit(port['normal_source']);F=shoulder_pads.pad_basis(n)
        rows=[]
        for angle in range(0,180,5):
            radians=np.deg2rad(angle);u=F[:,0]*np.cos(radians)+F[:,1]*np.sin(radians)
            holes=[]
            for sign in [-1,1]:
                centre=I+u*(sign*7.)
                observations={}
                for key,bone in bones.items():
                    observations[key]=[dict(point_source_mm=list(hit.position),
                        normal_mm=float((np.array(hit.position)-I)@n)) for hit in bone.ray_cast(centre+n*10.,centre-n*80.)]
                crossings=observations[target]
                holes.append(dict(sign=sign,axis_source_mm=centre.tolist(),bones=observations,
                    target_crossings=len(crossings),stock_mm=(crossings[0]['normal_mm']-crossings[1]['normal_mm']) if len(crossings)==2 else None))
            feasible=all(h['target_crossings']==2 for h in holes)
            rows.append(dict(angle_deg=angle,holes=holes,both_axes_cross_target_once=feasible,
                minimum_stock_mm=min(h['stock_mm'] for h in holes) if feasible else None,
                rear_surface_difference_mm=abs(holes[0]['bones'][target][1]['normal_mm']-holes[1]['bones'][target][1]['normal_mm']) if feasible else None))
        groups.append(dict(head=head,target=target,cap_inner_source_mm=I.tolist(),normal_source=n.tolist(),rows=rows))
        print(head,'eligible angles',[row['angle_deg'] for row in rows if row['both_axes_cross_target_once']],flush=True)
    save(DEST/'hole_axis_study.json',dict(scope=__doc__,parent_report=file.relative_to(ROOT).as_posix(),parent_sha256=sha(file),
        parent_verification_sha256=sha(check),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),
        screw_pitch_mm=14.,angle_step_deg=5.,groups=groups,no_mount_geometry_generated=True,stage_pass=False))


if __name__=='__main__':main()
