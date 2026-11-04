"""Source-faithful full-muscle orientations and machine-independent slice probes.

0.4mm nozzle / 0.2mm layers are analysis assumptions, NOT paper Table S2
or a selected manufacturing process. All heater temperatures are disabled.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'mujoco_arm/tools'))
import mesh_geom

DEST = ROOT/'prototype/s1_infraspinatus_trial/slicing_probe'


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    DEST.mkdir(exist_ok=True)
    trial = read(DEST.parent/'manifest.json')
    muscle = next(row for row in trial['parts'] if row['key'] == 'tpu:infraspinatus')
    source = ROOT/muscle['source_file']
    if sha(source) != muscle['source_sha256']:
        raise ValueError('Changed full-muscle source')
    selection = read(ROOT/'design/S1_CURRENT_CANDIDATE.json')
    graph_file = ROOT/selection['Y_pair_generation_report']
    if sha(graph_file) != selection['Y_pair_generation_sha256']:
        raise ValueError('Changed source frame definition')
    graph = next(row for row in read(graph_file)['attempts'] if row.get('selected') and row['module'] == 'infraspinatus')
    stations = graph['transported_frames']['stations']
    frame = stations[len(stations)//2]
    basis = np.array(frame['paper_basis_columns_source'])
    centre = np.array(frame['centre_source_mm'])
    vertices, faces = mesh_geom.read_stl(source)
    # Do not round/weld source coordinates. Reindex exact shared points,
    # keeping the same triangle order and double precision pose arithmetic.
    vertices, inverse = np.unique(vertices, axis=0, return_inverse=True)
    faces = inverse[faces]
    world_vertices = vertices+np.array(muscle['source_to_assembly_translation_mm'])
    rotations = {
        'paper_z_up': np.eye(3),
        'paper_z_down': np.diag([-1., 1., -1.]),
        'paper_y_up': np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]]),
        'paper_x_up': np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]]),
    }
    orientations = []
    for name, rotation in rotations.items():
        R = rotation@basis.T
        assert np.allclose(R.T@R, np.eye(3)) and np.linalg.det(R) > .99999
        transformed = (world_vertices-centre)@R.T
        lower, upper = transformed.min(axis=0), transformed.max(axis=0)
        shift = np.array([-(lower[0]+upper[0])/2, -(lower[1]+upper[1])/2, -lower[2]])
        positioned = transformed+shift
        folder = DEST/name
        folder.mkdir(exist_ok=True)
        target = folder/'full_muscle.obj'
        with target.open('w', encoding='ascii', newline='\n') as stream:
            stream.write('# Full current muscle, source triangles preserved, double precision orientation probe\n')
            for point in positioned:
                stream.write('v '+ ' '.join(format(float(value), '.17g') for value in point)+'\n')
            for face in faces:
                stream.write('f '+ ' '.join(str(int(value)+1) for value in face)+'\n')
        restored, read_faces = mesh_geom.load_obj(str(target))
        assert np.array_equal(read_faces, faces)
        assert len(read_faces) == 81668
        error = float(np.max(np.abs(restored-positioned)))
        triangles = restored[read_faces]
        cross = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
        area = np.linalg.norm(cross, axis=1)/2
        normals = cross/np.maximum(np.linalg.norm(cross, axis=1)[:, None], 1e-30)
        bottom = np.max(triangles[:, :, 2], axis=1) < .01
        down = (normals[:, 2] < -np.cos(np.pi/4)) & ~bottom
        orientations.append(dict(name=name, file=target.relative_to(ROOT).as_posix(), sha256=sha(target),
            triangles=len(faces), dimensions_mm=(upper-lower).tolist(),
            source_to_print_affine=np.c_[R, shift-R@centre].tolist(),
            print_to_source_affine=np.c_[R.T, centre-R.T@shift].tolist(),
            maximum_export_coordinate_error_mm=error,
            downward_surface_area_beyond_45deg_mm2=float(area[down].sum()),
            overhang_bin_excludes_bed=True,
            overhang_bin_is_not_support_volume_or_printability=True))
    common = '''# ANALYSIS ONLY. Not an executable manufacturing prescription.
# Assumptions: not paper Table S2, not selected hardware/material.
printer_technology = FFF
nozzle_diameter = 0.4
filament_diameter = 1.75
layer_height = 0.2
first_layer_height = 0.2
perimeters = 3
perimeter_generator = arachne
fill_density = 100%
fill_pattern = rectilinear
top_solid_layers = 4
bottom_solid_layers = 4
temperature = 0
first_layer_temperature = 0
bed_temperature = 0
first_layer_bed_temperature = 0
start_gcode = ; ANALYSIS ONLY - NOT FOR PRINTING
end_gcode = ; ANALYSIS ONLY - NOT FOR PRINTING
before_layer_gcode = ; ANALYSIS_LAYER_{layer_num} Z={layer_z}
bed_shape = 0x0,360x0,360x360,0x360
max_print_height = 360
gcode_flavor = reprap
use_relative_e_distances = 1
gcode_comments = 1
binary_gcode = 0
complete_objects = 0
skirts = 0
brim_width = 0
extrusion_width = 0.45
first_layer_extrusion_width = 0.45
external_perimeter_extrusion_width = 0.45
perimeter_extrusion_width = 0.45
infill_extrusion_width = 0.45
solid_infill_extrusion_width = 0.45
perimeter_speed = 20
external_perimeter_speed = 20
infill_speed = 20
solid_infill_speed = 20
top_solid_infill_speed = 20
bridge_speed = 15
support_material_speed = 20
first_layer_speed = 15
travel_speed = 120
support_material_auto = 1
support_material_threshold = 45
support_material_style = snug
support_material_contact_distance = 0.2
support_material_xy_spacing = 0.35
support_material_interface_layers = 2
support_material_buildplate_only = 0
cooling = 0
fan_always_on = 0
disable_fan_first_layers = 3
'''
    for name, support in [('no_support_analysis.ini', 0), ('support_analysis.ini', 1)]:
        (DEST/name).write_text(common+f'support_material = {support}\n', encoding='utf8', newline='\n')
    no_interface = common.replace('support_material_interface_layers = 2','support_material_interface_layers = 0')
    (DEST/'support_no_interface_analysis.ini').write_text(no_interface+'support_material = 1\nsupport_material_bottom_interface_layers = 0\nsupport_material_interface_contact_loops = 0\n',encoding='utf8',newline='\n')
    manifest = dict(scope=__doc__, muscle_source=muscle['source_file'], muscle_sha256=muscle['source_sha256'],
        frame_source=graph_file.relative_to(ROOT).as_posix(), frame_source_sha256=sha(graph_file),
        frame_station_index=len(stations)//2, paper_axes='centre station frame; muscle transported frames vary along its length',
        original_candidate_unchanged=True, orientations=orientations,
        parameters_source='explicit geometry analysis assumptions only; printer/material still unselected',
        not_verified=['actual nozzle/material, slicing profile or paper Table S2',
                      'supports removable, steel channel threading, strength, deformation or actuation'],
        stage_pass=False, actual_drive_range_deg=None)
    (DEST/'orientation_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf8', newline='\n')
    print(json.dumps([{k: row[k] for k in ['name', 'dimensions_mm', 'downward_surface_area_beyond_45deg_mm2']}
                      for row in orientations], indent=2))


if __name__ == '__main__':
    main()
