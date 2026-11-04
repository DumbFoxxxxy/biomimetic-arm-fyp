"""Find continuous rigid anterior-cap approaches without changing sources.

All current source objects remain present, including existing wire references.
The six complete planned routes are checked separately: a route conflict means
the cap must precede that wire, whose actual feeding is still unverified. This
checks a bounded corridor only. Never skip a blocked interval.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component,shoulder_pads

DEST=A.DESIGN/'s1_anterior_mount'


def main():
    file=DEST/'assembly_check.json';r=read(file)
    proof_file=DEST/'iterations/normal_cap_entry/verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(file) and proof['static_fit_and_contact_pass']
    solids={k:load_component(s) for k,s in r['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=read(route_file)['routes']
    route=next(x for x in routes if x['head']=='FJ1468');port=route['ports']['origin']
    n=A.D.unit(port['normal_source']);I=np.array(port['inner_face_source_mm']);F=shoulder_pads.pad_basis(n)
    key='cap:FJ1468/origin';moving=solids[key]
    # Screws and washers are fitted later, exactly as in the staged assembly.
    later={k for k in solids if k.startswith('origin:FJ1468/') and k!='origin:FJ1468/front_support'}
    fixed=set(solids)-later-{key};rows=[]
    for elevation in [0.,15.,30.,45.,60.,75.]:
        for azimuth in range(0,360,15):
            phi=np.deg2rad(azimuth);theta=np.deg2rad(elevation)
            direction=(F[:,0]*np.cos(phi)+F[:,1]*np.sin(phi))*np.cos(theta)+n*np.sin(theta)
            hit=None;count=0
            for distance in np.arange(0.,12.0001,.5):
                posed=moving.translate(direction*distance);count+=1
                for target in sorted(fixed):
                    v=float(A.I.intersection_volume(posed,solids[target]))
                    if v>.001:
                        hit=dict(offset_mm=float(distance),target=target,volume_mm3=v);break
                if hit:break
            rows.append(dict(azimuth_deg=azimuth,elevation_deg=elevation,direction_source=direction.tolist(),
                poses_checked=count,first_collision=hit,clear_to_12mm=hit is None))
        print('elevation',elevation,'clear',sum(x['clear_to_12mm'] for x in rows if x['elevation_deg']==elevation),flush=True)
    eligible=[x for x in rows if x['clear_to_12mm']]
    selected=max(eligible,key=lambda x:x['elevation_deg']) if eligible else None
    fine=[];wire_hits=[]
    if selected:
        direction=np.array(selected['direction_source'])
        for way,distances in [('in',np.arange(12.,-.0001,-.1)),('out',np.arange(0.,12.0001,.1))]:
            for distance in distances:
                posed=moving.translate(direction*distance)
                for target in sorted(fixed):
                    v=float(A.I.intersection_volume(posed,solids[target]))
                    if v>.001:fine.append(dict(way=way,offset_mm=float(distance),target=target,volume_mm3=v))
        for route in routes:
            corridor=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3)
            samples=[]
            for distance in np.arange(0.,12.0001,.1):
                v=float(A.I.intersection_volume(moving.translate(direction*distance),corridor))
                if v>.001:samples.append(dict(offset_mm=float(distance),volume_mm3=v))
            if samples:wire_hits.append(dict(head=route['head'],samples=samples))
    save(DEST/'cap_entry_study.json',dict(scope=__doc__,assembly_sha256=sha(file),
        verification_report=proof_file.relative_to(ROOT).as_posix(),verification_sha256=sha(proof_file),
        script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        fixed_keys=sorted(fixed),later_installed_keys=sorted(later),angle_step_deg=15.,coarse_step_mm=.5,
        corridor_length_mm=12.,rows=rows,selected=selected,selected_fine_step_mm=.1,
        selected_fine_roundtrip_poses=242 if selected else 0,selected_fine_hits=fine,
        planned_wire_corridor_hits=wire_hits,cap_must_precede_wires=[x['head'] for x in wire_hits],
        physical_wire_feeding_verified=False,
        local_cap_entry_pass=bool(selected) and not fine,whole_assembly_verified=False,stage_pass=False))
    print('selected',selected,'fine hits',len(fine),flush=True)


if __name__=='__main__':main()
