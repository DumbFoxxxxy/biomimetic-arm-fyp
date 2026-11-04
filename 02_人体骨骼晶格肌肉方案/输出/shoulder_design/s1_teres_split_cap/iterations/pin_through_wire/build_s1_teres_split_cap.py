"""Split distal teres cap to close sideways around the fixed muscle outlet.

Integral locating pins align halves, one M3 per half. No strength, loaded
pin retention, manufacture or whole-module installation claim.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,shoulder_pads,sha,read,save,export,load_component

DEST=A.DESIGN/'s1_teres_split_cap'


def main():
    DEST.mkdir(exist_ok=True)
    parent_file=A.DESIGN/'s1_teres_module_fit/assembly_check.json';parent=read(parent_file)
    assert parent['no_new_or_increased_overlap']
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    route=next(row for row in read(route_file)['routes'] if row['head']=='FJ1507')
    I=np.array(route['ports']['distal']['inner_face_source_mm']);n=A.D.unit(route['ports']['distal']['normal_source'])
    vector=np.array(route['core_path_source_mm'][-1])-I
    tangent=A.D.unit(vector-n*(vector@n));normal=A.D.unit(np.cross(tangent,n))
    F=np.stack([normal,tangent,n],axis=1);mf=A.I._manifold()
    cap_key='cap:FJ1507/distal';cap=load_component(parent['sources'][cap_key])
    local=cap.transform(np.c_[F.T,-F.T@I])
    left=local^mf.Manifold.cube([19.9,40.,20.]).translate([-20.,-20.,-5.])
    right=local^mf.Manifold.cube([19.9,40.,20.]).translate([.1,-20.,-5.])
    for y in [-6.5,6.5]:
        pin=rod([-1.1,y,2.],[1.2,y,2.],.75)
        bore=rod([-1.4,y,2.],[1.4,y,2.],.85)
        left=left+pin;right=right-bore
    parts={'cap:FJ1507/distal_A':left.transform(np.c_[F,I]),
           'cap:FJ1507/distal_B':right.transform(np.c_[F,I])}
    sources={k:v for k,v in parent['sources'].items() if k!=cap_key}
    solids={key:load_component(row) for key,row in sources.items()}
    for key,solid in parts.items():
        solids[key],sources[key]=export(solid,DEST/(key.split('/')[-1]+'.obj'),I)
    hits=[];count=0;keys=list(solids)
    for index,a in enumerate(keys):
        for b in keys[index+1:]:
            if a not in parts and b not in parts:continue
            count+=1;volume=float(A.I.intersection_volume(solids[a],solids[b]))
            if volume>.001:hits.append(dict(a=a,b=b,volume_mm3=volume))
    united=solids['cap:FJ1507/distal_A']+solids['cap:FJ1507/distal_B']
    report=dict(scope=__doc__,sources=sources,n_instances=len(sources),parent_report=parent_file.relative_to(ROOT).as_posix(),
        parent_sha256=sha(parent_file),script_sha256=sha(Path(__file__)),routes_sha256=sha(route_file),
        retired_keys=[cap_key],new_keys=sorted(parts),all_other_206_sources_unchanged=True,
        local_origin_source_mm=I.tolist(),local_basis_columns=F.tolist(),
        dimensions=dict(split_gap_mm=.2,pin_diameter_mm=1.5,pocket_diameter_mm=1.7,
            pin_tangent_offsets_mm=[-6.5,6.5],pin_axial_z_mm=2.,pin_normal_interval_mm=[-1.1,1.2],
            pocket_normal_interval_mm=[-1.4,1.4]),
        changed_pairs=count,static_hits=hits,static_geometry_pass=not hits,
        cap_material_added_mm3=float((united-cap).volume()),cap_material_removed_mm3=float((cap-united).volume()),
        free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        intended_body_owners={**{k:v for k,v in parent['intended_body_owners'].items() if k!=cap_key},
            **{k:'humerus_ideal_fastened' for k in parts}},
        adopted_as_current=False,whole_module_installation_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'assembly_check.json',report)
    print(dict(parts=len(solids),changed_pairs=count,hits=hits,cap_added=report['cap_material_added_mm3']),flush=True)
    if hits:raise ValueError('Split cap conflicts with source assembly; saved candidate is not acceptance')


if __name__=='__main__':main()
