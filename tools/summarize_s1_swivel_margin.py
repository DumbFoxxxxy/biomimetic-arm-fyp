"""Summarize the torque required to displace each cable-only swivel equilibrium."""
from build_s1_teres_origin_mount import ROOT, read, save, sha


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    source = folder/'extensor_late_reaction_audit.json'
    joint = next(row['combined_key_and_seat'] for row in read(source)['planar_bearing_trials']
                 if row.get('combined_key_and_seat', {}).get('sources'))
    data = joint['revised_elbow_motion']['tracked_swivels']['rows']
    results = {}
    for muscle, rows in data.items():
        groups = {}
        for row in rows:
            for sample in row['local_cable_only_equilibrium'].get('sensitivity', []):
                groups.setdefault(sample['offset_deg'], []).append({**sample, 'elbow_deg': row['elbow_deg']})
        limits = []
        for offset, samples in sorted(groups.items()):
            weakest = min(samples, key=lambda sample:
                          sample['balancing_moment_per_newton_Nmm'] * (1 if offset > 0 else -1))
            limits.append(dict(offset_deg=offset, poses=len(samples), limiting_elbow_deg=weakest['elbow_deg'],
                               opposing_moment_per_newton_Nmm=weakest['balancing_moment_per_newton_Nmm'],
                               peak_own_rotor_wire_overlap_mm3=max(item['wire_overlap_mm3'] for item in samples)))
        results[muscle] = limits
    save(folder/'swivel_resistance_margin.json', dict(
        source_sha256=sha(source), results=results,
        scope='Discrete local torque and own-rotor cable geometry. Scale moments by actual taut cable tension. Other contacts, friction hysteresis and TPU stiffness unverified.'))
    for muscle, limits in results.items():
        print(muscle, [(row['offset_deg'], round(row['opposing_moment_per_newton_Nmm'], 4),
                        row['limiting_elbow_deg']) for row in limits])


if __name__ == '__main__':
    main()
