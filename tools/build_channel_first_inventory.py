"""Record existing route evidence and unresolved thorax-to-elbow drive connections."""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def provenance(path):
    return dict(file=path.relative_to(ROOT).as_posix(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def path_metrics(points):
    points = np.asarray(points, dtype=float)
    assert points.ndim == 2 and points.shape[1] == 3 and np.isfinite(points).all()
    return dict(points_source_mm=points.tolist(), point_count=len(points),
                nominal_polyline_length_mm=float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum()))


def main():
    design = ROOT/'02_人体骨骼晶格肌肉方案/输出'
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    shoulder_file = ROOT/selection['assembly_report']
    assert provenance(shoulder_file)['sha256'] == selection['assembly_sha256']
    shoulder = read(shoulder_file)
    routes_file = design/'shoulder_design/s1_port_path_revision/routes_compact_distal_candidate.json'
    strand_file = design/'shoulder_design/s1_y2_transition_stock/free_strand_routes/route_audit.json'
    elbow_file = design/'integration_review/s1_elbow_shared_humerus/anatomical_elbow_baseline.json'
    baseline_assembly_file = elbow_file.parent/'assembly_check.json'
    elbow = read(elbow_file)
    assert elbow['original_assembly_sha256'] == provenance(baseline_assembly_file)['sha256']
    strands = read(strand_file)
    assert strands['assembly_sha256'] == provenance(shoulder_file)['sha256']
    rows = []
    owner = {'FJ1468': 'clavicle', 'FJ1513': 'scapula', 'FJ1467': 'scapula',
             'FJ1500': 'scapula', 'FJ1504': 'scapula', 'FJ1507': 'scapula'}
    for route in read(routes_file)['routes']:
        head = route['head']
        row = dict(id=head, region='shoulder', module=route['module'],
                   anatomical_origin_body=owner[head], distal_body='humerus',
                   local_route_reference=provenance(routes_file),
                   local_route=path_metrics(route['full_wire_path_source_mm']),
                   endpoint_ports_reference=route['ports'],
                   motor_mount=None, thorax_supply_channel=None,
                   proximal_conduit_reaction='existing local origin port candidate; full supply segment unresolved',
                   guide_support='TPU core channel and rigid end transitions; moving support law unresolved',
                   dynamic_route_validated=False)
        if head in ('FJ1500', 'FJ1504'):
            strand = next(item for item in strands['heads'] if item['head'] == head)
            neutral = strand['neutral']['route']
            row['local_route_status'] = 'Historical path reference; Y2 free terminal changes route topology'
            row['current_free_strand'] = dict(
                source=provenance(strand_file), endpoint_A_source_mm=neutral['exit_A_source_mm'],
                endpoint_B_source_mm=neutral['exit_B_source_mm'],
                neutral_chord_mm=neutral['exit_chord_mm'],
                installed_length_measured=neutral['installed_length_measured'],
                support='Free tension span between terminal exits; curved numerical family needs physical justification')
        else:
            row['local_route_status'] = 'Inherited candidate requiring verification against new whole-arm motion'
        rows.append(row)
    for label, muscle in elbow['muscles'].items():
        definition_file = ROOT/muscle['definition']
        assert provenance(definition_file)['sha256'] == muscle['definition_sha256']
        rows.append(dict(
            id=muscle['anatomical_id'], region='elbow', module=label,
            anatomical_origin_body='humerus', distal_body='forearm',
            local_route_reference=provenance(elbow_file),
            local_route=path_metrics(muscle['wire_path_source_mm']),
            pad_centres_source_mm=muscle['pad_centres_source_mm'],
            source_definition=provenance(definition_file), motor_mount=None, thorax_supply_channel=None,
            proximal_conduit_reaction='proximal rigid cap conduit seat candidate',
            distal_termination='distal cap crimp candidate; grip and retention unverified',
            guide_support='continuous TPU wall candidate connected to both end pads; compliance unresolved',
            dynamic_route_validated=False))
    report = dict(
        status='input_inventory_not_selected_global_route_design', coordinate_frame='BodyParts3D source mm',
        shoulder_selection=provenance(selection_file), shoulder_assembly=provenance(shoulder_file),
        original_elbow_assembly=provenance(baseline_assembly_file),
        shoulder_owners=shoulder['intended_body_owners'], routes=rows,
        missing_global_segments=[
            'thorax motor and conduit-reaction mounting locations',
            'supply routes from thorax to each shoulder and elbow local origin',
            'moving support law and reaction path for every TPU channel',
            'physical SC/AC or scapula-to-thorax mobility policy beyond fixed-context S1',
            'one selected whole-arm skeleton and collision inventory',
            'unilateral cable actuation and passive TPU spring/damper dynamic model'],
        final_support_policy=dict(
            grounded_body='thorax attached to tabletop base',
            motor_and_channel_reaction_mounts='thorax or explicitly assigned skeletal body',
            historical_test_fixture='engineering backboard and its clamps',
            tabletop_base_and_thorax_attachment_verified=False),
        historical_optional_interfaces=['cam and belt tracking', 'rotating pickup', 'deltoid contact relief'],
        complete=False)
    DEST.mkdir(exist_ok=True)
    (DEST/'route_inventory.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    lines = ['# 通道优先设计输入清单', '',
             '当前先完成胸廓—S1—肘部的作用路径、钢丝通道及反力支承，再设计完整 TPU 肌腹。',
             '本清单整理既有证据，八条局部路径尚未组成已选定的全局驱动方案。', '',
             '最终基准为胸廓通过底座固定在桌面。驱动器及通道反力支承以胸廓或指定骨段为基准。',
             '现有工程背板和夹具属于历史台架条件，基于它们的供线候选仅保留追溯。', '',
             '| 编号 | 区域 | 起始骨段 | 终止骨段 | 参考路径长度 / mm |',
             '| --- | --- | --- | --- | ---: |']
    for row in rows:
        lines.append(f"| {row['id']} | {row['region']} | {row['anatomical_origin_body']} | "
                     f"{row['distal_body']} | {row['local_route']['nominal_polyline_length_mm']:.2f} |")
    lines += ['', '长度沿各自源文件的路径计算，包含的端部引线和自由段不同，不能直接用作电机行程。',
              'FJ1500 和 FJ1504 另记录 S1-Y2 自由肌腱两端坐标，旧完整路径保留为参考。',
              '', '胸廓驱动器布置、供线段、通道随动及反力支承需要优先补齐。',
              '当前 S1 含肩胛／锁骨固定条件，胸廓到肩胛的主动运动还需要单独设计。',
              '凸轮、转座取角和肌腹让位暂列候选，采用与否由新的等效动力学结果决定。',
              '', '[坐标、支承状态和源文件指纹](route_inventory.json)', '']
    (DEST/'README.md').write_text('\n'.join(lines), encoding='utf-8')
    print('Inventoried', len(rows), 'local routes; no global dynamic route selected.')


if __name__ == '__main__':
    main()
