"""Report cycle evidence without treating free-space motion as acceptance."""
import csv
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance


def main():
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/elbow_extension_stop'
    check=read(folder/'model_check.json');parent=read(ROOT/check['parent']['file'])
    endpoints=read(ROOT/parent['endpoint_report']['file']);ids=[c['id'] for c in endpoints['channels']]
    rows=[];sources=[]
    labels=['k1e6','k1e7','k1e8']
    if (folder/'cycle_buffer_k1e6.json').exists():labels.append('buffer_k1e6')
    for label in labels:
        path=folder/f'cycle_{label}.json';r=read(path);sources.append(provenance(path))
        audit_path=folder/f'cycle_{label}_solid_screen.json'
        audit=read(audit_path) if audit_path.exists() else None
        if audit:assert audit['response']==provenance(path)
        channels=[]
        for head in ('FJ1486','FJ1477'):
            index=ids.index(head);trace=r['trace']
            lengths=np.array([s['tendon_lengths_mm'][index] for s in trace])
            max_turn=0.;min_radius=float('inf')
            for state in trace:
                points=np.array(state['wire_nodes_source_mm'][head]);segments=np.diff(points,axis=0)
                norms=np.linalg.norm(segments,axis=1);unit=segments/norms[:,None]
                turns=np.arccos(np.clip(np.sum(unit[:-1]*unit[1:],axis=1),-1.,1.))
                max_turn=max(max_turn,float(np.degrees(turns).max()))
                cross=np.linalg.norm(np.cross(segments[:-1],segments[1:]),axis=1)
                denominator=2*cross;valid=denominator>1e-12
                chord=np.linalg.norm(points[2:]-points[:-2],axis=1)
                radius=norms[:-1][valid]*norms[1:][valid]*chord[valid]/denominator[valid]
                if len(radius):min_radius=min(min_radius,float(radius.min()))
            channels.append(dict(id=head,peak_tension_N=r['peak_tensions_N'][index],
                centreline_length_range_mm=[float(lengths.min()),float(lengths.max())],
                centreline_excursion_mm=float(np.ptp(lengths)),slack_saved_states=sum(s['slack'][index] for s in trace),
                maximum_discrete_node_turn_deg=max_turn,minimum_three_node_radius_mm=min_radius))
        rows.append(dict(label=label,contact_solref=r['runtime_contact_parameters']['pair_solref'],
            contact_margin_mm=[m*1000 for m in r['runtime_contact_parameters'].get('pair_margin_m',[0.])],
            elbow_range_deg=r['elbow_range_deg'],final_elbow_deg=r['final_elbow_deg'],
            command_stroke_mm=1000*max(r['command']['pull_m']),
            maximum_contact_penetration_mm=r['maximum_contact_penetration_m']*1000,
            maximum_endpoint_constraint_error_mm=r['maximum_endpoint_constraint_error_m']*1000,
            saved_states=len(r['trace']),pad_edge_reached=r['pad_edge_reached'],pad_bottomed=r['pad_bottomed'],
            channels=channels,solid_audit=provenance(audit_path) if audit else None,
            first_rigid_interference=audit['first_rigid_interference'] if audit else 'pending',
            wire_interferences=audit['first_interferences'] if audit else 'pending',
            audit_rigid_states_checked=audit['rigid_states_checked'] if audit else 0,
            statistics_scope='Raw numerical whole-cycle values; motion interpretation ends at first solid interference.'))
    save(folder/'cycle_summary.json',dict(inputs=sources,trials=rows,
        generator=provenance(ROOT/'tools/summarize_channel_first_elbow_cycles.py'),
        scope='2s command pull/hold/payout/settle; not full K range. Frictionless cable and stop. Discrete bend radius from three beam nodes is a geometry warning metric, not continuous steel curvature or bend qualification. Missing audits remain pending; any first solid interference stops reachable-motion interpretation. No physical performance claim.'))
    with (folder/'cycle_summary.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['trial','elbow_min_deg','elbow_max_deg','final_elbow_deg','motor_pull_mm',
            'stop_penetration_mm','muscle','peak_tension_N','path_excursion_mm','min_discrete_radius_mm','rigid_audit_states','first_rigid_interference'])
        for r in rows:
            for c in r['channels']:writer.writerow([r['label'],*r['elbow_range_deg'],r['final_elbow_deg'],r['command_stroke_mm'],
                r['maximum_contact_penetration_mm'],c['id'],c['peak_tension_N'],c['centreline_excursion_mm'],
                c['minimum_three_node_radius_mm'],r['audit_rigid_states_checked'],str(r['first_rigid_interference'])])
    for row in rows:print(row['label'],row['channels'],'rigid audit',row['first_rigid_interference'],flush=True)


if __name__=='__main__':main()
