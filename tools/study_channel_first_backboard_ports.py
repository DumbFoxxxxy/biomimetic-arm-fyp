"""Screen explicit backboard pass-through openings for blocked supply paths."""
from itertools import combinations

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save
from build_channel_first_inventory import provenance
from screen_channel_first_distributed_supply import sample_curve


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    distributed_file = folder/'supply_channel_distributed_screen.json'
    report = read(distributed_file)
    assembly_file = ROOT/report['assembly']['file']
    assert provenance(assembly_file) == report['assembly']
    sources = read(assembly_file)['sources']
    board = load_component(sources['bench:backboard'])
    obstacles = {key: load_component(source) for key, source in sources.items()
                 if not key.startswith('tpu:') and key != 'bench:backboard'}
    fixed = {key: solid for key, solid in obstacles.items() if key.startswith('bench:')}
    accepted = {}
    for head, row in report['selected'].items():
        path, _ = sample_curve(np.asarray(row['controls_source_mm']), 801)
        accepted[head] = A.D.outer_envelope(path, 1.6)
    selected = {}
    bores = []
    for head in ('FJ1500', 'FJ1513'):
        candidates = [row for row in report['trials'][head]
                      if set(row.get('obstacle_hits_mm3', {})) == {'bench:backboard'}]
        candidates.sort(key=lambda row: row['obstacle_hits_mm3']['bench:backboard'])
        for candidate in candidates:
            path, bend_radius = sample_curve(np.asarray(candidate['controls_source_mm']), 801)
            if bend_radius < 8.:
                continue
            channel = A.D.outer_envelope(path, 1.6)
            if min((A.I.gap(channel, part, 2.) for part in accepted.values()), default=2.) < .5:
                continue
            if any(A.I.intersection_volume(channel, part) > .001 for part in obstacles.values()):
                continue
            bore = A.D.outer_envelope(path, 1.9)
            removal = board ^ bore
            hits = {key: volume for key, part in fixed.items()
                    if (volume := A.I.intersection_volume(bore, part)) > .001}
            if hits:
                continue
            bores.append(bore)
            accepted[head] = channel
            selected[head] = dict(controls_source_mm=candidate['controls_source_mm'],
                                  board_stock_removed_mm3=float(removal.volume()),
                                  refined_bend_radius_mm=bend_radius,
                                  required_support='lined backboard passage and retained conduit; final fitting pending')
            break
    mf = A.I._manifold()
    result = dict(distributed_screen=provenance(distributed_file), selected=selected,
                  bore_radius_mm=1.9, conduit_outer_radius_mm=1.6, sampled_points=801,
                  scope='Neutral backboard passage candidate. Machining, insertion, conduit restraint, '
                        'moving routes and board strength remain unverified.')
    if bores:
        revised = (board-mf.Manifold.batch_boolean(bores, mf.OpType.Add)).as_original().simplify(1e-6)
        try:
            actual, result['board_source'] = export(
                revised, folder/'backboard_supply_passage_candidate.obj',
                np.asarray(board.bounding_box()).reshape(2, 3).mean(axis=0))
        except (AssertionError, ValueError) as error:
            result['export_error'] = str(error) or type(error).__name__
        else:
            result['board_stock_removed_mm3'] = float((board-actual).volume())
            result['channel_to_revised_board_overlap_mm3'] = {
                head: A.I.intersection_volume(part, actual) for head, part in accepted.items()}
            result['supply_pair_checks'] = [
                dict(left=left, right=right, gap_mm_capped2=A.I.gap(accepted[left], accepted[right], 2.),
                     overlap_mm3=A.I.intersection_volume(accepted[left], accepted[right]))
                for left, right in combinations(accepted, 2)]
    save(folder/'backboard_supply_passage_screen.json', result)
    print('additional paths', list(selected), 'removed', result.get('board_stock_removed_mm3'),
          'export error', result.get('export_error'),
          'board overlaps', result.get('channel_to_revised_board_overlap_mm3'))


if __name__ == '__main__':
    main()
