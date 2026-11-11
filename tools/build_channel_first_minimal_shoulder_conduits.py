"""Generate hollow route prototypes, distinct from full muscles and approved mounts."""
from itertools import combinations

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save
from build_channel_first_inventory import provenance
from channel_first_tube_geometry import smooth_route, hollow_sweep


def simplify_path(path, tolerance=.01):
    """RDP within a stated centreline deviation; retain endpoint locations."""
    if len(path) <= 2:
        return path
    direction = path[-1]-path[0]
    length2 = direction@direction
    t = np.clip((path-path[0])@direction/max(length2, 1e-20), 0., 1.)
    distance = np.linalg.norm(path-(path[0]+t[:, None]*direction), axis=1)
    index = int(np.argmax(distance))
    if distance[index] <= tolerance:
        return path[[0, -1]]
    return np.vstack([simplify_path(path[:index+1], tolerance)[:-1],
                      simplify_path(path[index:], tolerance)])


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    audit_file = folder/'curved_shoulder_neutral_clearance.json'
    reference_file = folder/'thorax_arm_reference.json'
    reference = read(reference_file)
    context = dict(reference['sources'])
    context.update(read(folder/'endpoint_neutral_clearance.json')['cervical_context'])
    anatomy = {key: load_component(source) for key, source in context.items()}
    objects = {}
    rows = []
    rejected = []
    generated_solids, generated_wires = {}, {}
    for row in read(audit_file)['routes']:
        if not row['current_topology_compatible']:
            continue
        sequence = []
        for name in ('origin_transition', 'core', 'distal_transition'):
            for point in row['segments'][name]['points_source_mm']:
                if not sequence or np.linalg.norm(np.array(point)-sequence[-1]) > 1e-6:
                    sequence.append(point)
        original = np.asarray(sequence)
        for iterations in (144, 256, 400, 576):
            path, tangents, smoothing = smooth_route(original, iterations)
            if smoothing['minimum_sampled_bend_radius_mm'] >= 3.2:
                break
        else:
            rejected.append(dict(id=row['id'], reason='Provisional sampled bend radius below 3.2 mm'))
            continue
        conduit, origin = hollow_sweep(path, tangents)
        destination = folder/'minimal_shoulder_conduits'/f"{row['id']}.obj"
        destination.parent.mkdir(exist_ok=True)
        try:
            _, source = export(conduit, destination, np.zeros(3))
        except (AssertionError, ValueError) as error:
            rejected.append(dict(id=row['id'], reason=str(error) or
                'Exact coordinate welding creates nonmanifold edges; hollow tube not approved',
                status='geometry_generation_requires_revision'))
            destination.unlink(missing_ok=True)
            continue
        source['translation_source_mm'] = origin.tolist()
        # Check actual exported material in assembled source coordinates, not the old route.
        actual = load_component(source)
        wire_path = np.vstack([path[0]-tangents[0]*2., simplify_path(path, .005),
                               path[-1]+tangents[-1]*2.])
        wire = A.D.outer_envelope(wire_path, .3)
        passage_overlap = A.I.intersection_volume(wire, actual)
        assert passage_overlap < .001, 'Nominal cable intersects actual channel material'
        hits, wire_hits = {}, {}
        gap, nearest = 2., None
        for key, bone in anatomy.items():
            volume = A.I.intersection_volume(actual, bone)
            if volume > .001:
                hits[key] = volume
            wire_volume = A.I.intersection_volume(wire, bone)
            if wire_volume > .001:
                wire_hits[key] = wire_volume
            value = A.I.gap(actual, bone, gap)
            if value < gap:
                gap, nearest = value, key
        total_turn = float(np.arccos(np.clip(np.sum(tangents[:-1]*tangents[1:], axis=1), -1., 1.)).sum())
        objects[f"channel:{row['id']}"] = source
        generated_solids[row['id']], generated_wires[row['id']] = actual, wire
        rows.append(dict(id=row['id'], source=source, path_source_mm=path.tolist(),
            original_path_source_mm=original.tolist(), smoothing=smoothing,
            endpoint_tangents_source=tangents[[0, -1]].tolist(),
            outer_radius_mm=1.6, lumen_radius_mm=.8, nominal_wall_mm=.8,
            cross_section_sides=48, centreline_station_step_parameter_mm=.5,
            wire_diameter_reference_mm=.6, wire_access_extension_each_end_mm=2.,
            wire_channel_material_overlap_mm3=passage_overlap,
            channel_bone_overlap_mm3=hits, extended_wire_bone_overlap_mm3=wire_hits,
            actual_channel_gap_mm_capped2=gap, nearest_part=nearest,
            checked_anatomy_parts=len(anatomy),
            sampled_total_turn_radians=total_turn,
            anatomy_neutral_screen_clear=not hits and not wire_hits,
            length_mm=float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum()),
            end_opening_centres_source_mm=[path[0].tolist(), path[-1].tolist()],
            status='hollow_geometry_candidate_not_print_or_assembly_approved',
            support='Continuous TPU tube between inherited end-transition locations; end mounting, reaction support and deforming route law unresolved',
            contact_review='Actual changed sweep screened against anatomy; endpoint fittings, movement and soft contact still pending'))
        print(row['id'], 'bend radius', round(smoothing['minimum_sampled_bend_radius_mm'], 3),
              'tube/bone', list(hits), 'wire/bone', list(wire_hits), flush=True)
    pairs = []
    for left, right in combinations(generated_solids, 2):
        pairs.append(dict(left=left, right=right,
            channel_material_overlap_mm3=A.I.intersection_volume(generated_solids[left], generated_solids[right]),
            left_wire_right_channel_overlap_mm3=A.I.intersection_volume(generated_wires[left], generated_solids[right]),
            right_wire_left_channel_overlap_mm3=A.I.intersection_volume(generated_wires[right], generated_solids[left]),
            channel_gap_mm_capped2=A.I.gap(generated_solids[left], generated_solids[right], 2.)))
    save(folder/'minimal_shoulder_conduits.json', dict(
        route_audit=provenance(audit_file), anatomy_reference=provenance(reference_file),
        cervical_context_source=provenance(folder/'endpoint_neutral_clearance.json'),
        generation_script=provenance(ROOT/'tools/build_channel_first_minimal_shoulder_conduits.py'),
        sweep_script=provenance(ROOT/'tools/channel_first_tube_geometry.py'),
        conduits=rows, rejected_geometry=rejected, sources=objects, conduit_pair_screen=pairs,
        assumptions='Dimensions and 3.2 mm bend-radius screen are uncalibrated engineering candidates. Smooth indexed annular sweep, 48 sides, nominal centreline station parameter step <=0.5 mm. No complete muscle belly or lattice. Openings and nominal 0.6 mm wire passage checked geometrically, not physically threaded. Curvature is sampled at <=0.1 mm parameter spacing, not a mathematical global guarantee.',
        superseded_construction='Historical capsule subtraction had coordinate-weld failures in three guides. Smoothed indexed sweeps change the route and tangent directions, so original interface fit and motion approvals do not transfer.',
        exclusions=['Y2-modified FJ1500 and FJ1504', 'motor supply', 'mounts', 'force validation', 'deformation validation']))
    reference['sources'].update(objects)
    reference['status'] = 'thorax_arm_and_smoothed_hollow_channel_candidates_static_only'
    save(folder/'minimal_shoulder_conduits_scene.json', reference)
    print('Exported', len(rows), 'hollow conduits with exact OBJ round-trip checks;',
          len(rejected), 'rejected and withheld')


if __name__ == '__main__':
    main()
