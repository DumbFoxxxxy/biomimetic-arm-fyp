"""Replay actual extrusion segments and compare supports with the steel path.

These are nominal toolpath envelopes, not printed bead geometry or a proof
of removability. Supports touching the nominal steel core must be reviewed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'prototype/s1_infraspinatus_trial/slicing_probe'


def read_segments(path):
    position = np.zeros(3)
    absolute_e = True
    extrusion = 0.
    role, width, height = 'unclassified', .45, .2
    segments = []
    heating = []
    with path.open(encoding='utf8') as stream:
        for raw in stream:
            for prefix, name in [(';TYPE:', 'role'), (';WIDTH:', 'width'), (';HEIGHT:', 'height')]:
                if raw.startswith(prefix):
                    value = raw.strip()[len(prefix):]
                    if name == 'role': role = value
                    elif name == 'width': width = float(value)
                    else: height = float(value)
            fields = raw.split(';', 1)[0].split()
            if not fields: continue
            command = fields[0]
            if command in ['M104', 'M109', 'M140', 'M190']:
                values = [float(token[1:]) for token in fields[1:] if token.startswith('S')]
                heating.extend(value for value in values if value != 0)
            if command == 'G91': raise ValueError('Relative XYZ requires a separate explicit replay')
            if command == 'M82': absolute_e = True
            if command == 'M83': absolute_e = False
            if command == 'G92':
                for token in fields[1:]:
                    if token.startswith('E'): extrusion = float(token[1:])
            if command not in ['G0', 'G1']: continue
            values = {token[0]: float(token[1:]) for token in fields[1:] if token[0] in 'XYZEF'}
            target = np.array([values.get(axis, position[index]) for index, axis in enumerate('XYZ')])
            amount = 0.
            if 'E' in values:
                amount = values['E']-extrusion if absolute_e else values['E']
                extrusion = values['E'] if absolute_e else extrusion+values['E']
            if amount > 0 and np.linalg.norm(target[:2]-position[:2]) > 1e-9:
                segments.append((position.copy(), target.copy(), role, width, height, amount))
            position = target
    assert not heating, 'Analysis-only profile must have no nonzero heating commands'
    return segments


def distances_to_path(starts, ends, path):
    """Exact segment-distance candidates, with parallel cases covered by ends."""
    answer = np.full(len(starts), np.inf)
    u = ends-starts
    a = np.einsum('ij,ij->i', u, u)
    for p, q in zip(path[:-1], path[1:]):
        v = q-p
        c = float(v@v)
        if c < 1e-20: continue
        w = starts-p
        b = u@v
        d = np.einsum('ij,ij->i', u, w)
        e = w@v
        denominator = a*c-b*b
        valid = denominator > 1e-14
        s, t = np.zeros(len(a)), np.zeros(len(a))
        s[valid] = (b[valid]*e[valid]-c*d[valid])/denominator[valid]
        t[valid] = (a[valid]*e[valid]-b[valid]*d[valid])/denominator[valid]
        inside = valid & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
        candidate = w+s[:, None]*u-t[:, None]*v
        distance = np.full(len(a), np.inf)
        distance[inside] = np.linalg.norm(candidate[inside], axis=1)
        for point in [starts, ends]:
            projection = np.clip(((point-p)@v)/c, 0., 1.)
            distance = np.minimum(distance, np.linalg.norm(point-(p+projection[:, None]*v), axis=1))
        for point in [p, q]:
            projection = np.clip(np.einsum('ij,ij->i', point-starts, u)/a, 0., 1.)
            distance = np.minimum(distance, np.linalg.norm(point-(starts+projection[:, None]*u), axis=1))
        answer = np.minimum(answer, distance)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orientation', required=True)
    parser.add_argument('--support', action='store_true')
    parser.add_argument('--block-channel', action='store_true')
    parser.add_argument('--blocker-radius',type=float,default=1.25)
    parser.add_argument('--no-interface',action='store_true')
    parser.add_argument('--vertical-shadow',action='store_true')
    args = parser.parse_args()
    kind = 'with_channel_blocker' if args.block_channel else 'with_support' if args.support else 'no_support'
    if args.block_channel and args.blocker_radius!=1.25:kind+='_r'+str(args.blocker_radius).replace('.','p')
    if args.vertical_shadow:kind+='_shadow'
    if args.no_interface:kind+='_no_interface'
    folder = DEST/args.orientation/kind
    report = json.loads((folder/'slice_report.json').read_text(encoding='utf8'))
    gcode = ROOT/report['raw_gcode_file']
    assert hashlib.sha256(gcode.read_bytes()).hexdigest() == report['raw_gcode_sha256']
    segments = read_segments(gcode)
    total_volume = sum(row[5]*math.pi*(1.75/2)**2 for row in segments)
    summary_volume = sum(row['volume_mm3'] for row in report['gcode_summary']['roles'].values())
    assert abs(total_volume-summary_volume) < 1e-5, 'Independent extrusion replay disagrees with summary'
    manifest = json.loads((DEST/'orientation_manifest.json').read_text(encoding='utf8'))
    pose = next(row for row in manifest['orientations'] if row['name'] == args.orientation)
    trial = json.loads((DEST.parent/'manifest.json').read_text(encoding='utf8'))
    route = trial['active_route_source_mm']
    full = np.array(route['full_wire_path_source_mm'])
    i0 = int(np.argmin(np.linalg.norm(full-np.array(route['ports']['origin']['inner_face_source_mm']), axis=1)))
    i1 = int(np.argmin(np.linalg.norm(full-np.array(route['ports']['distal']['inner_face_source_mm']), axis=1)))
    assert i0 < i1
    transform = np.array(pose['source_to_print_affine'])
    path = full[i0:i1+1]@transform[:, :3].T+transform[:, 3]+np.array(report['slicer_instance_translation_mm'])
    supports = [row for row in segments if row[2].lower().startswith('support')]
    support_result = dict(segment_count=len(supports), minimum_centerline_distance_mm=None,
                          centerlines_inside_nominal_wire_radius_count=0, possible_bead_envelope_overlap_count=0)
    selected_z = []
    if supports:
        starts, ends = np.array([row[0] for row in supports]), np.array([row[1] for row in supports])
        widths, heights = np.array([row[3] for row in supports]), np.array([row[4] for row in supports])
        starts[:, 2] -= heights/2
        ends[:, 2] -= heights/2
        distances = distances_to_path(starts, ends, path)
        inside = distances < .3
        # Radius encloses the declared width/height rectangle. This is an
        # overestimate of nominal bead extent, with no sag/material model.
        radius = np.hypot(widths/2, heights/2)
        possible = distances < .3+radius
        idx = np.flatnonzero(inside)
        support_result.update(minimum_centerline_distance_mm=float(distances.min()),
            centerlines_inside_nominal_wire_radius_count=int(inside.sum()),
            possible_bead_envelope_overlap_count=int(possible.sum()),
            interpretation='0.3mm nominal steel radius; conservative rectangular-bead enclosing capsule; not physical obstruction/removability proof',
            example_segments=[dict(start_mm=starts[i].tolist(), end_mm=ends[i].tolist(),
                distance_mm=float(distances[i]), width_mm=float(widths[i]), height_mm=float(heights[i]), role=supports[i][2])
                for i in idx[:12]])
        selected_z = [round(float(supports[i][1][2]), 6) for i in idx]
    zvalues = sorted({round(float(row[1][2]), 6) for row in segments})
    critical = max(set(selected_z), key=selected_z.count) if selected_z else zvalues[len(zvalues)//2]
    reference = support_result
    baseline_file = DEST/args.orientation/'with_support/channel_path_check.json'
    if args.block_channel and baseline_file.exists():
        reference = json.loads(baseline_file.read_text(encoding='utf8'))['support_to_steel_core']
    examples = reference.get('example_segments', [])
    anchor = next((row for row in examples if row['height_mm'] < .3), examples[0] if examples else None)
    roi = None
    if anchor:
        requested = anchor['start_mm'][2]+anchor['height_mm']/2
        critical = min(zvalues, key=lambda value: abs(value-requested))
        roi = (np.array(anchor['start_mm'][:2])+np.array(anchor['end_mm'][:2]))/2
    levels = list(dict.fromkeys([zvalues[0], critical, zvalues[len(zvalues)//2], zvalues[3*len(zvalues)//4]]))
    # Native vector toolpath preview, drawn from actual positive extrusion
    # records. No generated appearance is used as evidence of deformation.
    colours = {'support': '#7e8b9d', 'interface': '#356ca1', 'body': '#bd4937'}
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="850" viewBox="0 0 1120 850"><rect width="1120" height="850" fill="white"/>',
           '<g font-family="Arial,sans-serif" fill="#233047"><text x="28" y="34" font-size="22">Actual PrusaSlicer extrusion paths — '+args.orientation+' / '+kind+'</text>',
           '<text x="28" y="62" font-size="14">Red: body | grey: support | blue: support interface | black circles: nominal steel-core section approximation</text>',
           '<text x="28" y="84" font-size="14">Analysis assumptions: 0.4mm nozzle, 0.2mm layers. Profiles/material unverified. Not for printing.</text></g>']
    dense = []
    for p, q in zip(path[:-1], path[1:]):
        dense.extend(p+(q-p)*np.linspace(0., 1., max(2, int(np.linalg.norm(q-p)/.1)+1))[:, None])
    dense = np.array(dense)
    for index, level in enumerate(levels[:4]):
        rows = [row for row in segments if abs(row[1][2]-level) < 1e-5]
        pts = np.vstack([np.array([row[0] for row in rows]), np.array([row[1] for row in rows])])
        lower, upper = pts[:, :2].min(axis=0)-1., pts[:, :2].max(axis=0)+1.
        if index == 1 and roi is not None:
            lower, upper = roi-6., roi+6.
        scale = min(490/(upper[0]-lower[0]), 290/(upper[1]-lower[1]))
        cx, cy = 28+(index%2)*550, 112+(index//2)*358
        def xy(p): return cx+18+(p[0]-lower[0])*scale, cy+320-(p[1]-lower[1])*scale
        label=' | 12mm bore-region zoom' if index==1 and roi is not None else ''
        svg.append(f'<rect x="{cx}" y="{cy}" width="525" height="342" fill="#fcfcfa" stroke="#d5d9dd"/><text x="{cx+12}" y="{cy+22}" font-family="Arial" font-size="14">Z={level:.2f} mm{label} | {len(rows)} segments</text>')
        svg.append(f'<defs><clipPath id="panel{index}"><rect x="{cx+12}" y="{cy+32}" width="501" height="299"/></clipPath></defs><g clip-path="url(#panel{index})">')
        for row in rows:
            x1, y1 = xy(row[0]); x2, y2 = xy(row[1])
            colour = colours['interface'] if 'support' in row[2].lower() and 'interface' in row[2].lower() else colours['support'] if row[2].lower().startswith('support') else colours['body']
            svg.append(f'<path d="M{x1:.3f},{y1:.3f}L{x2:.3f},{y2:.3f}" stroke="{colour}" stroke-width="{max(.3,row[3]*scale):.3f}" stroke-linecap="round"/>')
        for p in dense[np.abs(dense[:, 2]-(level-.1)) < .3][::2]:
            x, y = xy(p)
            if cx < x < cx+525 and cy+30 < y < cy+342:
                r = math.sqrt(max(0., .3**2-(p[2]-(level-.1))**2))*scale
                svg.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" fill="none" stroke="black" stroke-width=".5"/>')
        svg.append('</g>')
    svg.append('</svg>')
    (folder/'toolpath_layers.svg').write_text(''.join(svg), encoding='utf8')
    result = dict(scope=__doc__, slice_report_sha256=hashlib.sha256((folder/'slice_report.json').read_bytes()).hexdigest(),
        independent_positive_extrusion_segments=len(segments), total_positive_extrusion_volume_mm3=total_volume,
        agrees_with_job_summary=True, all_temperatures_disabled=True,
        internal_channel_source_point_range=[i0, i1], support_to_steel_core=support_result,
        preview_Z_levels_mm=levels, preview_zoom_anchor_from_automatic_supports=anchor,
        analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        toolpath_preview_sha256=hashlib.sha256((folder/'toolpath_layers.svg').read_bytes()).hexdigest(),
        support_removal_verified=False, physical_print_verified=False,
        stage_pass=False, actual_drive_range_deg=None)
    (folder/'channel_path_check.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(args.orientation, kind, json.dumps(support_result), flush=True)


if __name__ == '__main__':
    main()
