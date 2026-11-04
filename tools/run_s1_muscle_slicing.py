"""Run actual PrusaSlicer offline probes and record positive extrusion roles.

Raw G-code is stored in ignored tmp only, with heating disabled. It is NOT
a manufacturing handoff. Printer, material and paper Table S2 are unknown.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'prototype/s1_infraspinatus_trial/slicing_probe'
DEFAULT_SLICER = ROOT/'.deps/prusaslicer-2.9.6/PrusaSlicer-2.9.6/prusa-slicer-console.exe'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gcode_summary(path):
    xyz = [0., 0., 0.]
    e = 0.
    absolute_e = True
    role = 'unclassified'
    layers = set()
    roles = defaultdict(lambda: dict(segments=0, xy_length_mm=0., filament_mm=0., volume_mm3=0.))
    minimum_width, maximum_width = float('inf'), 0.
    width = .45
    commands = defaultdict(int)
    with path.open(encoding='utf8', errors='replace') as stream:
        for raw in stream:
            if raw.startswith(';TYPE:'):
                role = raw.strip().split(':', 1)[1]
            if raw.startswith(';WIDTH:'):
                width = float(raw.strip().split(':', 1)[1])
            code = raw.split(';', 1)[0].strip()
            if not code:
                continue
            parts = code.split()
            command = parts[0]
            commands[command] += 1
            if command == 'M82':
                absolute_e = True
            elif command == 'M83':
                absolute_e = False
            elif command == 'G92':
                for token in parts[1:]:
                    if token.startswith('E'):
                        e = float(token[1:])
            elif command in ['G2', 'G3']:
                raise ValueError('Arc toolpaths need explicit parser support; no silent linear approximation')
            elif command in ['G0', 'G1']:
                values = {letter: float(value) for letter, value in re.findall(r'([XYZEF])([-+]?(?:\d+(?:\.\d*)?|\.\d+))', code)}
                new = [values.get(letter, xyz[index]) for index, letter in enumerate('XYZ')]
                amount = 0.
                if 'E' in values:
                    amount = values['E']-e if absolute_e else values['E']
                    e = values['E'] if absolute_e else e+values['E']
                length = math.hypot(new[0]-xyz[0], new[1]-xyz[1])
                if amount > 0 and length > 1e-9:
                    layers.add(round(new[2], 6))
                    row = roles[role]
                    row['segments'] += 1
                    row['xy_length_mm'] += length
                    row['filament_mm'] += amount
                    row['volume_mm3'] += amount*math.pi*(1.75/2)**2
                    minimum_width, maximum_width = min(minimum_width, width), max(maximum_width, width)
                xyz = new
    return dict(extrusion_Z_count=len(layers), min_extrusion_Z_mm=min(layers) if layers else None,
                max_extrusion_Z_mm=max(layers) if layers else None, roles=dict(roles),
                minimum_reported_width_mm=minimum_width if layers else None, maximum_reported_width_mm=maximum_width,
                command_counts=dict(commands),
                volume_note='Positive XY extrusion only, 1.75mm assumed filament; not measured material, mass or mechanical performance')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orientation', choices=['paper_z_up', 'paper_z_down', 'paper_y_up', 'paper_x_up'], required=True)
    parser.add_argument('--support', action='store_true')
    parser.add_argument('--block-channel', action='store_true')
    parser.add_argument('--blocker-radius',type=float,default=1.25)
    parser.add_argument('--no-interface',action='store_true')
    parser.add_argument('--vertical-shadow',action='store_true')
    parser.add_argument('--slicer', type=Path, default=DEFAULT_SLICER)
    args = parser.parse_args()
    manifest = json.loads((DEST/'orientation_manifest.json').read_text(encoding='utf8'))
    source = next(row for row in manifest['orientations'] if row['name'] == args.orientation)
    model = ROOT/source['file']
    assert sha(model) == source['sha256']
    if args.block_channel:
        assert args.support
        tag='r'+str(args.blocker_radius).replace('.','p')
        suffix='_shadow' if args.vertical_shadow else ''
        filename='blocker_manifest.json' if args.blocker_radius==1.25 and not suffix else f'blocker_{tag}{suffix}_manifest.json'
        blocker = json.loads((DEST/filename).read_text(encoding='utf8'))
        chosen = next(row for row in blocker['packages'] if row['orientation'] == args.orientation)
        model = ROOT/chosen['file']
        assert sha(model) == chosen['sha256']
    kind = 'with_channel_blocker' if args.block_channel else 'with_support' if args.support else 'no_support'
    if args.block_channel and args.blocker_radius!=1.25:kind+='_'+tag
    if args.vertical_shadow:
        assert args.block_channel
        kind+='_shadow'
    if args.no_interface:
        assert args.support
        kind+='_no_interface'
    profile = DEST/('support_no_interface_analysis.ini' if args.no_interface else 'support_analysis.ini' if args.support else 'no_support_analysis.ini')
    output = DEST/args.orientation/kind
    output.mkdir(exist_ok=True)
    scratch = ROOT/'tmp/s1_slicing_probe'/args.orientation/kind
    scratch.mkdir(parents=True, exist_ok=True)
    datadir = scratch/'slicer_datadir'
    gcode = scratch/'ANALYSIS_ONLY_NOT_FOR_PRINTING.gcode'
    common = [str(args.slicer.resolve()), '--datadir', str(datadir), '--load', str(profile),
              '--threads', '2', '--dont-arrange', '--center', '180,180']
    info = subprocess.run(common+['--info', str(model)], capture_output=True, text=True, encoding='utf8', errors='replace', timeout=60)
    (output/'import_info.txt').write_text(info.stdout+info.stderr, encoding='utf8', newline='\n')
    if info.returncode:
        raise ValueError('Slicer import inspection failed')
    command = common+['--export-gcode', '--save', str(output/'effective_config.ini'), '--output', str(gcode), str(model)]
    started = time.monotonic()
    try:
        result = subprocess.run(command, capture_output=True, text=True, encoding='utf8', errors='replace', timeout=600)
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as error:
        code, stdout, stderr = None, str(error.stdout or ''), str(error.stderr or '')+'\n600s process timeout; child reaped by subprocess.run'
    elapsed = time.monotonic()-started
    (output/'slicer_stdout.txt').write_text(stdout, encoding='utf8', newline='\n')
    (output/'slicer_stderr.txt').write_text(stderr, encoding='utf8', newline='\n')
    success = code == 0 and gcode.exists() and gcode.stat().st_size > 0
    report = dict(scope=__doc__, orientation=args.orientation, support_enabled=args.support,
        slicer_path=str(args.slicer.relative_to(ROOT)) if args.slicer.is_relative_to(ROOT) else str(args.slicer),
        slicer_sha256=sha(args.slicer), version='PrusaSlicer 2.9.6 project-local portable',
        release_archive_url='https://github.com/prusa3d/PrusaSlicer/releases/download/version_2.9.6/PrusaSlicer-2.9.6.zip',
        release_archive_sha256='5aaf22e42f95accecfa122d23a835911f289ecc2ff606db3e83d637ddcc0a209',
        input_file=model.relative_to(ROOT).as_posix(), input_sha256=sha(model), source_pose=source['file'],
        channel_blocker_enabled=args.block_channel, channel_blocker_radius_mm=args.blocker_radius if args.block_channel else None,
        support_interface_disabled=args.no_interface,
        vertical_channel_shadow=args.vertical_shadow,
        profile=profile.relative_to(ROOT).as_posix(), profile_sha256=sha(profile),
        command=command, elapsed_s=elapsed, exit_code=code, analysis_slicing_completed=success,
        import_info=info.stdout+info.stderr, warnings=[line for line in (stdout+'\n'+stderr).splitlines() if any(s in line.lower() for s in ['warning', 'repair', 'empty', 'error'])],
        raw_gcode_file=gcode.relative_to(ROOT).as_posix(), raw_gcode_stored_in_ignored_tmp=True,
        raw_gcode_sha256=sha(gcode) if success else None,
        slicer_instance_translation_mm=[180., 180., 0.],
        gcode_summary=gcode_summary(gcode) if success else None,
        actual_printer_and_material_selected=False, support_removal_verified=False, printability_verified=False,
        actual_drive_range_deg=None, stage_pass=False)
    (output/'slice_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps({key: report[key] for key in ['orientation', 'support_enabled', 'exit_code', 'analysis_slicing_completed', 'elapsed_s', 'warnings']}, indent=2), flush=True)
    if not success:
        raise ValueError('Actual slicing did not complete; see preserved engine logs')


if __name__ == '__main__':
    main()
