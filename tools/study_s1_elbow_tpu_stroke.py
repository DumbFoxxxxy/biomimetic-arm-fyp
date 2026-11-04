"""Measure cable and TPU span demands for the relocated elbow candidate."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha


def arc_length(points):
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def screen_extensor_anchor(start, guide, distal, centre, axis):
    centre, axis = np.array(centre), np.array(axis)
    axial_point = centre+axis*((distal-centre)@axis)
    radial = distal-axial_point
    radius = float(np.linalg.norm(radial))
    trials = []
    for shift in (0., 5., 10., 15., 20., 25.):
        if shift >= radius:
            break
        anchor = distal-radial*(shift/radius)
        angles = np.arange(121)
        lengths, moment_arms = [], []
        for angle in angles:
            transform = A.I.rotation_about(centre, axis, int(angle))
            far = transform[:, :3]@anchor+transform[:, 3]
            lengths.append(arc_length(np.stack([start, guide, far])))
            direction = A.D.unit(guide-far)
            moment_arms.append(float(np.cross(far-centre, direction)@axis))
        trials.append(dict(radial_shift_mm=shift, anchor_source_mm=anchor.tolist(),
                           axis_radius_mm=radius-shift,
                           samples=[dict(elbow_deg=int(angle), stroke_mm=length-lengths[0],
                                         signed_moment_arm_mm=arm)
                                    for angle, length, arm in zip(angles, lengths, moment_arms, strict=True)],
                           minimum_extensor_moment_arm_mm=-max(moment_arms)))
    return trials


def audit_anchor_geometry(base, trials, start, guide, centre, axis):
    assembly = base/'assembly_tpu_readback_candidate.json'
    sources = read(assembly)['sources']
    forearm = load_component(sources['bone:forearm'])
    humerus = load_component(read(base/'elbow_relocated_bone_seats/bone_seat_check.json')['source'])
    elbow_file = (ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces'
                  /'mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
                  /'continuous_channel_candidate/assembly.json')
    owners = {f"elbow:{row['key']}": row['owner'] for row in read(elbow_file)['instances']}
    station_file = base/'elbow_proximal_rotated_station_screen.json'
    groups = read(station_file)['moving_groups']
    pair = read(base/'elbow_relocated_channel_cores/channel_core_check.json')['choice']['pair']
    frame = read(A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json')['frame_source']
    shaft = np.array(frame['basis_columns'])[:, 2]
    ownership_conflicts = []
    for muscle, head in [('flexor', 'FJ1486'), ('extensor', 'FJ1477')]:
        pose = A.I.rotation_about(frame['origin_mm'], shaft, pair[f'{muscle}_angle_deg'])
        pose[:, 3] += shaft*pair[f'{muscle}_shift_mm']
        for key in groups[head]:
            owner = owners.get(key, '')
            if owner != 'elbow' and not owner.endswith('_distal_pivot'):
                continue
            part = load_component(sources[key])
            point = np.array(part.bounding_box()).reshape(2, 3).mean(axis=0)
            moved = pose[:, :3]@point+pose[:, 3]
            ownership_conflicts.append(dict(key=key, original_owner=owner,
                                            source_to_relocated_centre_mm=float(np.linalg.norm(moved-point))))
    joints = {key: load_component(row) for key, row in sources.items()
              if key.startswith('elbow:elbow_')}
    for row in trials:
        anchor = np.array(row['anchor_source_mm'])
        probe_radius = .1
        probe = A.I._manifold().Manifold.sphere(probe_radius, 24).translate(anchor)
        embedded = bool(A.I.intersection_volume(probe, forearm) > .001)
        row['anchor_probe_overlaps_forearm'] = embedded
        row['anchor_forearm_gap_mm_capped20'] = A.I.gap(probe, forearm, 20.)
        hits = []
        for angle in range(121):
            pose = A.I.rotation_about(centre, axis, angle)
            end = pose[:, :3]@anchor+pose[:, 3]
            cable = A.D.outer_envelope(np.stack([start, guide, end]), .3)
            obstacles = {'bone:humerus': humerus, 'bone:forearm': forearm.transform(pose)}
            obstacles.update({key: part.transform(pose) if owners.get(key) == 'elbow' else part
                              for key, part in joints.items()})
            collisions = {key: volume for key, part in obstacles.items()
                          if (volume := float(A.I.intersection_volume(cable, part))) > .001}
            if collisions:
                hits.append(dict(elbow_deg=angle, hits=collisions))
        row['bone_and_joint_wire_hits'] = hits
    selected = next(row for row in trials if row['radial_shift_mm'] == 15.)
    anchor = np.array(selected['anchor_source_mm'])
    existing = load_component(sources['elbow:FJ1477_active_wire_crimp_barrel'])
    root = np.array(existing.bounding_box()).reshape(2, 3).mean(axis=0)
    feed = A.D.unit(anchor-guide)
    side = A.D.unit(np.cross(feed, anchor-root))
    support_trials = []
    mf = A.I._manifold()
    # In forearm coordinates the anchor is fixed and the incoming cable
    # direction sweeps. Only the local stock at the new terminal is relieved.
    fans = []
    for angle in range(121):
        inverse = A.I.rotation_about(centre, axis, -angle)
        local_guide = inverse[:, :3]@guide+inverse[:, 3]
        direction = A.D.unit(local_guide-anchor)
        fans.append(rod(anchor-direction*.4, anchor+direction*12., .55))
    clearance = mf.Manifold.batch_boolean(fans, mf.OpType.Add)
    for bow in (0., 3., -3., 6., -6.):
        mid = (root+anchor)/2+side*bow
        terminal = rod(anchor-feed*.15, anchor+feed*3.85, 1.5)
        arm = (rod(root, mid, 1.25)+rod(mid, anchor+feed, 1.25)
               +mf.Manifold.sphere(1.25, 24).translate(mid))
        sleeve = rod(anchor-feed*2, anchor+feed*6, .5)
        raw = (existing+arm+terminal)-sleeve
        candidate = (raw-clearance).as_original().simplify(1e-6)
        added = candidate-existing
        first = None
        for angle in range(121):
            pose = A.I.rotation_about(centre, axis, angle)
            end = pose[:, :3]@anchor+pose[:, 3]
            cable = A.D.outer_envelope(np.stack([start, guide, end]), .3)
            obstacles = {'bone:humerus': humerus, 'bone:forearm': forearm.transform(pose), 'wire': cable}
            obstacles.update({key: part.transform(pose) if owners.get(key) == 'elbow' else part
                              for key, part in joints.items()})
            moved = added.transform(pose)
            hits = {key: vol for key, part in obstacles.items()
                    if (vol := float(A.I.intersection_volume(moved, part))) > .001}
            if hits:
                first = dict(elbow_deg=angle, hits=hits,
                             terminal_wire_overlap_mm3=float(A.I.intersection_volume(
                                 cable, (terminal-sleeve).transform(pose))),
                             arm_wire_overlap_mm3=float(A.I.intersection_volume(
                                 cable, (arm-sleeve).transform(pose))))
                break
        trial = dict(bow_mm=bow, first_hit=first, components=len(candidate.decompose()),
                     swept_relief_removed_mm3=float((raw-candidate).volume()),
                     original_root_removed_mm3=float((existing-candidate).volume()),
                     relief_radius_mm=.55)
        support_trials.append(trial)
        if first is None and trial['components'] == 1:
            try:
                _, trial['source'] = export(candidate, base/'elbow_relocated_bone_seats/extensor_inboard_support.obj', anchor)
            except (AssertionError, ValueError) as error:
                trial['export_error'] = str(error)
                continue
            break
    return dict(assembly_sha256=sha(assembly), elbow_owner_sha256=sha(elbow_file),
                station_screen_sha256=sha(station_file), relocation_ownership_conflicts=ownership_conflicts,
                inboard_15mm_support_trials=support_trials,
                checked_joint_parts=sorted(joints), step_deg=1.,
                scope='Bone and central elbow hardware only. Relocated distal fittings and guide mounting still require design.')


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    bridge_file = base/'elbow_relocated_bridge_screen.json'
    core_file = base/'elbow_relocated_channel_cores/channel_core_check.json'
    bridge, core = read(bridge_file), read(core_file)
    waypoint = np.array(bridge['extensor_crimp_elbow_sweep']
                        ['extensor_fixed_guide_screen']['selected']['guide_source_mm'])
    centre, axis = A.frames.elbow_axis()
    rows = {}
    anchor_trials = []
    for channel in core['channels']:
        label = channel['label']
        start, end = np.array(core['choice'][f'{label}_endpoints_source_mm'])
        neutral = np.vstack([start, channel['path_source_mm'], end])
        station = np.r_[0., np.cumsum(np.linalg.norm(np.diff(neutral, axis=0), axis=1))]
        if label == 'extensor':
            route = np.stack([start, waypoint, end])
            anchor_trials = screen_extensor_anchor(start, waypoint, end, centre, axis)
            anchor_geometry = audit_anchor_geometry(base, anchor_trials, start, waypoint, centre, axis)
        else:
            route = neutral
        reference_wire = arc_length(route)
        # TPU sleeves stop 1 mm before the nominal anchors. On the extensor,
        # treat the fixed guide as the proximal end of the moving TPU span.
        proximal_trim = arc_length(route[:2]) if label == 'extensor' else 1.
        neutral_span = reference_wire-proximal_trim-1.
        samples = []
        for angle in range(121):
            transform = A.I.rotation_about(centre, axis, angle)
            far = transform[:, :3]@end+transform[:, 3]
            fractions = np.array([0., 0., 1.]) if label == 'extensor' else station/station[-1]
            posed = route+np.outer(fractions, far-end)
            wire_length = arc_length(posed)
            span = wire_length-proximal_trim-1.
            samples.append(dict(elbow_deg=angle, wire_length_mm=wire_length,
                                wire_stroke_mm=wire_length-reference_wire,
                                modelled_TPU_span_mm=span,
                                average_TPU_span_change_fraction=(span-neutral_span)/neutral_span))
        rows[label] = dict(neutral_wire_mm=reference_wire, neutral_TPU_span_mm=neutral_span,
                           assumed_distal_free_span_mm=1., assumed_proximal_inactive_span_mm=proximal_trim,
                           poses=samples)
    save(base/'elbow_tpu_stroke_demand.json', dict(
        bridge_sha256=sha(bridge_file), channel_sha256=sha(core_file), angle_step_deg=1.,
        rows=rows, extensor_anchor_trials=anchor_trials,
        anchor_geometry=anchor_geometry,
        scope='Kinematic centerline demands with fixed terminal offsets; no local TPU strain, load, slip or measured stroke.'))
    for label, row in rows.items():
        print(label, [row['poses'][angle] for angle in (30, 60, 90, 120)])


if __name__ == '__main__':
    main()
