"""Explicit stop-arm insertion and provisional nut-tool clearance probes."""
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,rod
from build_channel_first_inventory import provenance


def main():
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_extension_stop'
    assembly_file=folder/'assembly.json';assembly=read(assembly_file)
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    E=np.array(assembly['origin_source_mm']);F=np.array(assembly['frame_columns']);axis=F[:,0]
    sign=assembly['selected_side'];results=[]
    parent={k:v for k,v in solids.items() if not k.startswith('stop:')}
    # Existing joint is already assembled; arm/fasteners are installed in stages.
    stationary={**parent,**{k:v for k,v in solids.items() if k.startswith('stop:seat_washer')}}
    checked=0;first=None
    for distance in np.arange(30.,-.01,-1.):
        placed=solids['stop:arm'].translate(axis*sign*distance);hits=[]
        for key,solid in stationary.items():
            volume=float(A.I.intersection_volume(placed,solid))
            if volume>.001:hits.append(dict(target=key,overlap_mm3=volume))
        checked+=1
        if hits:first=dict(offset_mm=float(distance),hits=hits);break
    results.append(dict(stage='arm_axial_insertion_from_outboard',step_mm=1.,start_offset_mm=30.,
                        states_checked=checked,first_interference=first))
    # Tool dimensions are clearance hypotheses, not a selected commercial tool.
    # Probe a coaxial tubular 5.5AF nut tool (outer R4mm) from fork interior.
    for index,y in enumerate((15.,30.)):
        p=lambda x:E+F@np.array([sign*x,y,0.])
        sleeve=rod(p(12.),p(34.),4.)-rod(p(11.),p(35.),3.3)
        hits=[]
        for key,solid in solids.items():
            volume=float(A.I.intersection_volume(sleeve,solid))
            if volume>.001:hits.append(dict(target=key,overlap_mm3=volume))
        results.append(dict(stage=f'coaxial_nut_tool_{index}',outer_radius_mm=4.,inner_radius_mm=3.3,
                            axial_span_mm=[12.,34.],hits=hits))
    compact=[]
    for index,y in enumerate((15.,30.)):
        p=lambda x:E+F@np.array([sign*x,y,0.])
        for head_start,head_radius,shaft_radius in [(25.5,3.6,1.8),(25.5,4.,1.8),(25.5,4.5,1.8),
                                                    (26.2,3.6,1.8),(26.2,3.6,1.5)]:
            # Socket ends at rear washer, rather than passing through it.
            tool=(rod(p(head_start),p(33.5),head_radius)-rod(p(head_start-.2),p(34.),3.25))
            tool+=rod(p(12.),p(head_start),shaft_radius)
            hits=[]
            for key,solid in solids.items():
                volume=float(A.I.intersection_volume(tool,solid))
                if volume>.001:hits.append(dict(target=key,overlap_mm3=volume))
            compact.append(dict(index=index,y_mm=y,head_radius_mm=head_radius,
                                head_axial_mm=[head_start,33.5],shaft_radius_mm=shaft_radius,hits=hits))
    save(folder/'compact_tool_exploration.json',dict(assembly=provenance(assembly_file),
        generator=provenance(ROOT/'tools/audit_channel_first_elbow_stop_installation.py'),hypothetical_tools=compact,
        status='Stationary tool envelope exploration only; no real tool selected. Handle, insertion and hand access absent.'))
    save(folder/'installation_screen.json',dict(assembly=provenance(assembly_file),
        generator=provenance(ROOT/'tools/audit_channel_first_elbow_stop_installation.py'),stages=results,
        scope='Arm insertion: 1mm sampled path with bone/joint/cap solids and two preplaced seat washers, bolts/nuts not yet installed. Tool envelope is provisional and stationary only; it is not an insertion path or a confirmed real tool. Screwdriver/head socket geometry, washers/nuts insertion, main elbow insertion and human access remain unverified. No TPU envelope check.'))
    for row in results:print(row,flush=True)


if __name__=='__main__':main()
