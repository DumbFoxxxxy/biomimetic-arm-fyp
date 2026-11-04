"""Screen rotated posterior screw assemblies against TPU and active routes.

Cap and bone need planned bores and are not accepted by this preliminary
screen. The final builder must drill real geometry and check all pairs.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,M,rod,sha,read,save,load_component


def main():
    study_file=A.DESIGN/'s1_deltoid_origin_holes/hole_axis_study.json';study=read(study_file)
    parent=read(ROOT/study['parent_report']);assert sha(ROOT/study['parent_report'])==study['parent_sha256']
    group=next(row for row in study['groups'] if row['head']=='FJ1513');I=np.array(group['cap_inner_source_mm']);n=np.array(group['normal_source'])
    fixed={k:load_component(row) for k,row in parent['sources'].items() if k not in ['bone:scapula','cap:FJ1513/origin']}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    port=next(row for row in routes if row['head']=='FJ1513')['ports']['origin'];thickness=float((np.array(port['outer_face_source_mm'])-I)@n)
    wire_solids={row['head']:A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3) for row in routes}
    zone=rod(I-n*18.,I+n*.2,12.);rows=[]
    for row in group['rows']:
        if not row['both_axes_cross_target_once']:continue
        rear=min(h['bones']['bone:scapula'][1]['normal_mm'] for h in row['holes'])-2.2
        length=next(x for x in [20.,25.,30.,35.,40.] if x>=thickness+.75-(rear-3.15)+.5)
        pad=fixed['tpu:deltoid']
        for h in row['holes']:
            point=np.array(h['axis_source_mm']);pad=pad-(rod(point-n*45,point+n*10,1.7)^zone)
        hits=[]
        for i,h in enumerate(row['holes']):
            hardware=M.m3_parts(np.array(h['axis_source_mm']),n,thickness,rear,length)
            for name,solid in hardware.items():
                for key,obstacle in fixed.items():
                    volume=float(A.I.intersection_volume(solid,pad if key=='tpu:deltoid' else obstacle))
                    if volume>.001:hits.append(dict(index=i,part=name,target=key,volume_mm3=volume))
                for head,wire in wire_solids.items():
                    volume=float(A.I.intersection_volume(solid,wire))
                    if volume>.001:hits.append(dict(index=i,part=name,target='active_wire:'+head,volume_mm3=volume))
        rows.append(dict(angle_deg=row['angle_deg'],minimum_ray_stock_mm=row['minimum_stock_mm'],
            preliminary_fastener_hits=hits,preliminary_clear=not hits))
        print(row['angle_deg'],'stock',round(row['minimum_stock_mm'],3),'hits',len(hits),flush=True)
    save(study_file.parent/'posterior_fastener_study.json',dict(scope=__doc__,parent_sha256=study['parent_sha256'],
        axis_study_sha256=sha(study_file),routes_sha256=sha(route_file),script_sha256=sha(Path(__file__)),rows=rows,
        cap_and_bone_bore_fit_pending=True,backing_plate_pending=True,stage_pass=False))


if __name__=='__main__':main()
