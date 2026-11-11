"""Separate physical sheath seats from floating display-stub tips."""
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save,load_component,mesh_geom
from build_channel_first_inventory import provenance


def main():
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_supply_interface'
    assembly_file=folder/'assembly.json';assembly=read(assembly_file);rows=[]
    for r in assembly['interfaces']:
        head=r['id'];P=np.array(r['origin_source_mm']);axis=np.array(r['reaction_axis_source'])
        def top(source):
            vertices,_=mesh_geom.manifold_mesh_to_arrays(load_component(source).to_mesh64())
            return P+axis*float(np.max((vertices-P)@axis))
        outlet=top(assembly['sources'][f'supply:{head}_ferrule'])
        tip=top(assembly['conduit_stub_sources']['sheath:'+head])
        seat=np.array(r['sheath_bearing_point_source_mm'])
        rows.append(dict(id=head,owner='humerus',motor_mount=None,complete_supply_path=None,
            physical_ferrule_bearing_source_mm=r['ferrule_bearing_point_source_mm'],
            physical_sheath_reaction_seat_source_mm=seat.tolist(),
            physical_ferrule_outlet_source_mm=outlet.tolist(),outward_axis_source=axis.tolist(),
            ferrule_sheath_engagement_mm=float((outlet-seat)@axis),
            reference_only_stub_tip_source_mm=tip.tolist(),reference_tip_fixed=False,
            legacy_field_meaning='assembly.interfaces.supply_exit_source_mm is the neutral DISPLAY stub tip, not a rigid support or a fixed routing via point.',
            candidate_supply_boundary='Sheath compression seats at the actual ferrule inner step. Translation support exists there; a clamped tangent through4mm engagement is a future approximation requiring retention/contact verification. All outer sheath nodes must remain free or have explicit guide/contact supports.',
            steel_connection='Steel passes through sheath/ferrule/cap lumen into original muscle channel; proximal wire is not anchored to the cap. Remote spool and distal steel termination remain separate.',
            local_model_note='Internal beam root at original first guide is still clamped as an equivalent boundary, not a resolved lattice-to-cap joint or full external supply.',
            unresolved=['sheath-to-ferrule pullout retention','thread/preload/tool access','remote motor/reaction mount','supported full supply curve and bending/compression response','distal steel termination review']))
    save(folder/'route_endpoints.json',dict(assembly=provenance(assembly_file),
        generator=provenance(ROOT/'tools/describe_channel_first_elbow_ports.py'),
        coordinate_frame='BP3D source mm; original proportions and coordinates',routes=rows,
        scope='Authoritative roles for future supply routing on this local candidate; not selected whole-arm routes. No engineering backboard. Neutral display tip is explicitly not a support.'))
    for row in rows:print(row['id'],'real sheath seat',row['physical_sheath_reaction_seat_source_mm'],'engagement',row['ferrule_sheath_engagement_mm'],flush=True)


if __name__=='__main__':main()
