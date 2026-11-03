"""Separate two-platen y/z coupons, with and without an integrated wire tube.

Same 2x2x2 core and same print axes for each pair, no orthogonal grip frame.
The specimens are BCC-inspired, not author CAD or calibrated full muscles.
P0 is platen compression; P1 is a hand-wire assembly candidate with an actual
pressure washer and uncrimped stop. Neither has measured TPU performance.
"""
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import frames
import interference as I
import audit_shoulder_retained_joint_fit as A
from build_shoulder_unified_muscle import local_export
from build_shoulder_z_coupons import cell_pieces,rod

DEST=A.DESIGN/'s1_actuation_calibration/axis_coupons'
RECIPES={'reference_15_10_30':[15.,10.,30.],'candidate_10_8_10':[10.,8.,10.]}


def core(cell):
    pieces=[]
    for index in itertools.product(range(2),repeat=3):
        centre=(np.array(index)+.5)*cell
        pieces.extend(cell_pieces(centre,np.array(cell)/2))
    return I._manifold().Manifold.batch_boolean(pieces,I._manifold().OpType.Add)


def specimen(name,cell,axis,with_tube,body):
    folder=DEST/name;folder.mkdir(parents=True,exist_ok=True)
    index={'y':1,'z':2}[axis];direction=np.eye(3)[index]
    dims=np.array(cell)*2;length=float(dims[index]);centre=dims/2
    pieces=[body];thickness=3.
    for start in [-thickness,length]:
        low=np.full(3,-1.);size=dims+2
        low[index]=start;size[index]=thickness
        pieces.append(I._manifold().Manifold.cube(size).translate(low))
    tpu=I._manifold().Manifold.batch_boolean(pieces,I._manifold().OpType.Add)
    # Both pairs use x-separated M3 grip holes; there is no rigid/flexible
    # plate on a side parallel to the load axis to shunt the gauge deformation.
    for sign in [-1,1]:
        p=centre.copy();p[0]+=sign*min(8.,dims[0]/2-4.);p[index]=-4
        q=p.copy();q[index]=1
        tpu=tpu-rod(p,q,1.7)
    sources={};solids={}
    if with_tube:
        start=centre.copy();end=centre.copy();start[index]=-3;end[index]=length+3
        tube=rod(start,end,1.8)
        tpu=tpu+tube
        p=start-direction*2;q=end+direction*5
        tpu=tpu-rod(p,q,.86)
        cap_inner=end
        cap=(rod(cap_inner,cap_inner+direction*3,5.)-
             rod(cap_inner-direction,cap_inner+direction*4,.86)-
             rod(cap_inner+direction*.4,cap_inner+direction*4,1.4))
        stop=(rod(cap_inner+direction*.45,cap_inner+direction*2.95,1.25)-
              rod(cap_inner,cap_inner+direction*3.5,.45))
        wire=rod(start-direction*10,cap_inner+direction*3,.3)
        solids.update(pressure_washer=cap,uncrimped_stop=stop,nominal_wire_NOT_print=wire)
    solids['TPU']=tpu
    actual={}
    for key,solid in solids.items():
        actual[key],sources[key]=local_export(solid,folder/(key+'.stl'),np.zeros(3))
        if not sources[key]['mesh_checks_pass']:raise ValueError('Coupon mesh failed '+name+'/'+key)
    intersections=[]
    for a,b in itertools.combinations(actual,2):
        volume=I.intersection_volume(actual[a],actual[b])
        if volume>1e-9:intersections.append(dict(a=a,b=b,volume_mm3=volume))
    if any(row['volume_mm3']>.001 for row in intersections):raise ValueError('Nominal coupon material conflict '+name)
    stop_capture=None
    if with_tube:
        stop_capture={str(distance):I.intersection_volume(actual['uncrimped_stop'].translate(-direction*distance),actual['pressure_washer'])
                      for distance in [0.,.04,.06]}
    row=dict(name=name,cell_xyz_mm=list(cell),strut_diameter_mm=2.,cell_count_xyz=[2,2,2],
        paper_load_axis=axis,print_axes='unchanged XYZ for both load-axis variants; intended +Z up',
        gauge_length_mm=length,nominal_gauge_area_mm2=float(np.prod(dims)/length),
        grip_faces_only=[axis+'-',axis+'+'],orthogonal_grip_faces=False,
        continuous_channel=with_tube,channel_outer_radius_mm=1.8 if with_tube else None,
        channel_bore_radius_mm=.86 if with_tube else None,
        experiment='P1 hand-wire compression assembly, not a stiffness material test' if with_tube
                   else 'P0 two-platen compression; use nominal stress/strain for axis comparison',
        coordinate_frame='standalone coupon millimetres, NOT anatomical source coordinates',
        sources=sources,material_intersections=intersections,nominal_stop_capture=stop_capture,
        measured_shortening_mm=None,measured_force_N=None,calibrated_force_model=False,stage_pass=False)
    (folder/'design.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf8',newline='\n')
    print(name,'gauge',length,'area',row['nominal_gauge_area_mm2'],'parts',len(sources),'readback OK',flush=True)
    return row


def run():
    DEST.mkdir(parents=True,exist_ok=True);rows=[]
    for recipe,cell in RECIPES.items():
        body=core(np.array(cell))
        for axis,tube in itertools.product(['y','z'],[False,True]):
            rows.append(specimen(recipe+'_'+axis+('_tube' if tube else '_bare'),cell,axis,tube,body))
    report=dict(scope=__doc__,script_sha256=A.hashed(Path(__file__).resolve()),
        cell_geometry_script_sha256=A.hashed(Path(__file__).with_name('build_shoulder_z_coupons.py')),
        cell_is_author_CAD=False,material_candidate='TPU95A; actual printer, filament and slicing not selected',
        specimens_per_variant=3,geometry_variants=8,suggested_total_samples=24,coupons=rows,
        material_or_muscle_calibration_transferred=False,
        not_verified=['physical printing, supports, channel openness and effective load transfer',
            'compression, hysteresis, recovery and strain limits; no universal 25% assumption',
            'P1 cap is wire-preloaded pressure contact, not a tensile specimen or full anatomical muscle'],
        stage_pass=False)
    (DEST/'manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    # Never overwrite user measurements on regeneration.
    blank=DEST/'measurements_BLANK.json'
    if not blank.exists():
        blank.write_text(json.dumps(dict(scope='blank template; all readings are unmeasured',
            records=[dict(coupon=r['name'],specimen=i,printer=None,nozzle_mm=None,filament_batch=None,
                layer_mm=None,print_orientation=None,measured_gauge_mm=None,force_sensor=None,
                cycles=[],natural_length_before_mm=None,residual_60s_mm=None,notes=None)
                for r in rows for i in range(1,4)]),indent=2)+'\n',encoding='utf8',newline='\n')
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
