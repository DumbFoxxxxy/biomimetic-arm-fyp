"""Create a dimensioned stop drawing, detail scene and unmeasured grip record."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'prototype/s1_infraspinatus_trial/screw_stop_v1'


def main():
    design_file = DEST/'stop_design.json'
    design = json.loads(design_file.read_text(encoding='utf8'))
    dims = design['nominal_dimensions']
    assert dims['head_stock_diameter_mm'] == 4. and dims['head_ear_side_flat_local_x_mm'] == -1.5
    assert dims['head_z_interval_mm'] == [3.05, 7.05] and dims['stub_z_interval_mm'] == [.45, 3.05]
    assert dims['through_bore_mm'] == .8 and dims['radial_thread_axis_local_z_mm'] == 5.05
    digest = hashlib.sha256(design_file.read_bytes()).hexdigest()
    # Orthographic profile uses bottom-of-stub Z=0 (assembly local Z=.45).
    # Dimensions describe the nominal blank; the threaded port must be
    # machined to its callout, not drilled to the visualization envelope.
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="820" viewBox="0 0 1100 820">
<defs><marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#334155"/></marker>
<clipPath id="flat"><rect x="720" y="220" width="200" height="200"/></clipPath></defs>
<style>text{font-family:Arial,sans-serif;fill:#243347;font-size:16px}.dim{stroke:#334155;stroke-width:1;fill:none;marker-start:url(#arrow);marker-end:url(#arrow)}.ext{stroke:#8894a3;stroke-width:1}.body{fill:#f4e7bd;stroke:#27394f;stroke-width:2}.axis{stroke:#63748a;stroke-dasharray:8 5}.hidden{stroke:#3876a6;stroke-width:1.5;stroke-dasharray:5 4;fill:none}</style>
<rect width="1100" height="820" fill="white"/>
<text x="48" y="48" style="font-size:24px;font-weight:bold">FJ1500 reusable METAL cable stop</text>
<text x="48" y="78">Machining candidate | units: mm | nominal dimensions | no grip or load rating</text>
<text x="190" y="134">SIDE: local X-Z projection</text><text x="689" y="134">TOP: stock + ear-side flat</text>
<path class="body" d="M230 480L230 376L220 376L220 216L360 216L360 376L330 376L330 480Z"/>
<path class="axis" d="M280 175V510"/><path class="hidden" d="M264 216V480M296 216V480"/>
<path class="hidden" d="M280 254H400M280 338H400"/>
<path d="M292 278L304 256H412V336H304L292 314Z" fill="#cbd5e1" fill-opacity=".85" stroke="#465a72" stroke-width="1.5"/>
<path class="hidden" d="M380 278H412M380 314H412"/>
<path class="axis" d="M275 296H445"/>
<text x="425" y="270">M2 x 3 flat point</text><text x="425" y="294">0.9 AF hex socket</text><text x="425" y="318">nominal first touch</text>
<path class="ext" d="M220 216V175M360 216V175"/><path class="dim" d="M220 184H360"/><text x="271" y="174">3.5</text>
<path class="ext" d="M230 480H170M230 376H170"/><path class="dim" d="M180 480V376"/><text x="128" y="433">2.6</text>
<path class="ext" d="M220 216H170M220 376H170"/><path class="dim" d="M180 216V376"/><text x="128" y="270">4.0</text>
<path class="ext" d="M230 480V530M330 480V530"/><path class="dim" d="M230 518H330"/><text x="245" y="546">stub DIA 2.5</text>
<path class="dim" d="M110 480V296"/><path class="ext" d="M110 296H210M110 480H210"/><text x="62" y="386">4.6</text>
<text x="48" y="575">Axial hole DIA 0.8 THROUGH; bottom of stub is drawing Z=0.</text>
<circle cx="780" cy="320" r="80" fill="none" stroke="#a8b0bc" stroke-dasharray="5 5"/>
<circle class="body" cx="780" cy="320" r="80" clip-path="url(#flat)"/><path d="M720 267.085L720 372.915" stroke="#27394f" stroke-width="2"/>
<circle cx="780" cy="320" r="16" fill="white" stroke="#27394f" stroke-width="2"/>
<path class="axis" d="M680 320H915M780 220V420"/>
<path class="hidden" d="M780 278H860M780 362H860"/>
<path class="ext" d="M700 240V205M860 240V205"/><path class="dim" d="M700 214H860"/><text x="727" y="198">stock DIA 4.0</text>
<path class="dim" d="M700 450H720"/><path class="ext" d="M700 350V458M720 373V458"/><text x="660" y="480">flat depth 0.5 (X=-1.5)</text>
<text x="873" y="354">radial thread</text><text x="873" y="378">M2 x 0.4</text>
<path d="M720 320H650" class="ext"/><text x="617" y="296">ear side</text>
<text x="48" y="624">1. Machined METAL body; grade and manufacturing tolerances not selected. Do not print in TPU/PLA.</text>
<text x="48" y="654">2. Radial hole: M2 x 0.4 thread. The STL DIA 2.1 clearance is NOT a tap-drill instruction.</text>
<text x="48" y="684">3. Threads are not modelled. Confirm purchased screw tip/socket dimensions before manufacture.</text>
<text x="48" y="714">4. Assembly: retain original cap; nominal 0.05 axial seating gap; steel wire nominal DIA 0.6.</text>
<text x="48" y="744">5. First-touch geometry does not prove grip. Measure slip, force and strand damage on a separate sample.</text>
<text x="48" y="788" style="font-size:11px">Design SHA256: DESIGN_SHA</text>
</svg>'''.replace('DESIGN_SHA', digest)
    (DEST/'machining_candidate.svg').write_text(svg, encoding='utf8', newline='\n')
    full = json.loads((DEST/'full_shoulder_check.json').read_text(encoding='utf8'))
    keys = ['free_terminal:FJ1500', 'stop_screw:FJ1500/body', 'stop_screw:FJ1500/M2x3',
            'wire_active:FJ1500/extended_reference']
    (DEST/'detail_check.json').write_text(json.dumps(dict(
        scope='Actual cap, screw-stop and active steel; surrounding TPU/bone/tendons hidden; no force or deformation',
        sources={key: full['sources'][key] for key in keys}, n_instances=len(keys), stage_pass=False),
        ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    record = DEST/'grip_measurements_BLANK.json'
    if not record.exists():
        record.write_text(json.dumps(dict(scope='Unmeasured stop-only and then whole-muscle physical trial',
            design_sha256=digest, metal_grade=None, body_actual_dimensions_mm=None,
            actual_thread_fit=None, screw_supplier_and_part=None, screw_actual_tip_diameter_mm=None,
            wire_supplier_and_construction=None, wire_actual_diameter_mm=None, wire_coating=None,
            stop_axial_fixture_description=None, cap_seating_and_tool_access=None,
            tightening_method=None, tightening_torque_Nm=None, intended_trial_tension_N=None,
            cycles=[dict(cycle=index, pull_displacement_mm=None, peak_force_N=None,
                         wire_slip_mm=None, loosening=None, visible_wire_damage=None,
                         cap_damage=None, evidence_files=[]) for index in range(1, 6)],
            release_and_reclamp_success=None, axial_grip_verified=False, whole_muscle_pull_verified=False,
            actual_drive_range_deg=None, stage_pass=False), ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print('Machining SVG, four-object detail scene and unmeasured grip form prepared.')


if __name__ == '__main__':
    main()
