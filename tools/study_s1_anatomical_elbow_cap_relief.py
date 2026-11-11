"""Evaluate local cap relief while retaining the original elbow muscle geometry."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, rod, save, sha


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    assembly = folder/'assembly_check.json'
    baseline_file = folder/'anatomical_elbow_baseline.json'
    shoulder_file = folder/'anatomical_elbow_shoulder_clearance.json'
    baseline = read(baseline_file)
    shoulder = read(shoulder_file)
    assert baseline['original_assembly_sha256'] == sha(assembly)
    assert shoulder['baseline_sha256'] == sha(baseline_file)
    sources = read(assembly)['sources']
    revised = next(row['candidate_source'] for row in shoulder['trials'] if 'candidate_source' in row)
    deltoid = load_component(revised)
    prefix = 'elbow:FJ1477_proximal_'
    cap = load_component(sources[prefix+'rigid_cap'])
    extensor = load_component(baseline['muscles']['extensor']['collision_proxy'])
    adapter = load_component(sources['elbow:FJ1477_bone_adapter'])
    wire = A.D.outer_envelope(np.array(baseline['muscles']['extensor']['wire_path_source_mm']), 2.)
    mf = A.I._manifold()
    guards = {}
    for key, source in sources.items():
        if key.startswith(prefix) and ('washer' in key or 'bolt' in key):
            part = load_component(source)
            lo, hi = np.array(part.bounding_box()).reshape(2, 3)
            guards[key] = mf.Manifold.cube(hi-lo+.5).translate(lo-.25)
    rows = []
    relocation = []
    blind_trials = []
    ear_trials = []
    exported_interface = False
    for margin in (0., .05, .15):
        cuts = [cap ^ deltoid.translate(direction*margin)
                for direction in ([np.zeros(3)] if not margin else [np.zeros(3), *np.eye(3), *(-np.eye(3))])]
        candidate = (cap-mf.Manifold.batch_boolean(cuts, mf.OpType.Add)).as_original().simplify(1e-6)
        removed = cap-candidate
        hits = {key: volume for key, guard in guards.items()
                if (volume := float(A.I.intersection_volume(removed, guard))) > .001}
        row = dict(sample_offset_mm=margin, cap_removed_mm3=float(removed.volume()),
                   remaining_overlap_mm3=float(A.I.intersection_volume(candidate, deltoid)),
                   positive_components=sum(part.volume() > 1e-8 for part in candidate.decompose()),
                   fastener_guard_loss_mm3=hits,
                   muscle_interface_loss_mm3=float(A.I.intersection_volume(removed, extensor)),
                   wire_region_loss_mm3=float(A.I.intersection_volume(removed, wire)))
        rows.append(row)
        if margin == .05:
            group = {key: load_component(source) for key, source in sources.items()
                     if key.startswith(prefix) and key.endswith('_left')
                     and any(name in key for name in ('bolt', 'washer', 'nut'))}
            head = np.array(group[prefix+'head_washer_left'].bounding_box()).reshape(2, 3).mean(axis=0)
            tail = np.array(group[prefix+'nut_washer_left'].bounding_box()).reshape(2, 3).mean(axis=0)
            normal = A.D.unit(head-tail)
            u = A.D.unit(np.cross(normal, np.eye(3)[np.argmin(abs(normal))]))
            v = np.cross(normal, u)
            right = {key: load_component(source) for key, source in sources.items()
                     if key.startswith(prefix) and key.endswith('_right')}
            span = np.linalg.norm(head-tail)
            blockers = {'deltoid': deltoid, 'extensor': extensor, 'wire_guard': wire, **right}
            for radius in (2., 4., 6., 8., 10.):
                for degrees in range(0, 360, 30):
                    angle = np.radians(degrees)
                    shift = radius*(u*np.cos(angle)+v*np.sin(angle))
                    moved = {key: part.translate(shift) for key, part in group.items()}
                    hits = {name: peak for name, obstacle in blockers.items()
                            if (peak := max(float(A.I.intersection_volume(part, obstacle))
                                            for part in moved.values())) > .001}
                    core = rod(tail+shift, head+shift, 1.5)
                    thickness_proxy = float(A.I.intersection_volume(core, candidate))/(np.pi*1.5**2)
                    land = rod(tail+shift, head+shift, 3.5)-core
                    relocation.append(dict(radius_mm=radius, azimuth_deg=degrees, shift_source_mm=shift.tolist(),
                                           hardware_hits=hits, span_mm=float(span),
                                           cap_axis_stock_length_proxy_mm=thickness_proxy,
                                           cap_annular_stock_mm3=float(A.I.intersection_volume(land, candidate))))
            # Search the remaining cap stock with a shallow blind fastener
            # approached from the original head side at three inclinations.
            for shift in (4., 7., 10.):
                for degrees in range(0, 360, 45):
                    theta = np.radians(degrees)
                    radial = u*np.cos(theta)+v*np.sin(theta)
                    point = (head+tail)/2+radial*shift
                    for tilt in (30., 45., 60.):
                        alpha = np.radians(tilt)
                        approach = normal*np.cos(alpha)+radial*np.sin(alpha)
                        hits = candidate.ray_cast(point+approach*20., point-approach*20.)
                        if len(hits) < 2:
                            continue
                        face = np.array(hits[0].position)
                        depth = float(np.linalg.norm(face-np.array(hits[1].position)))
                        if depth < 3.5:
                            continue
                        screw = rod(face-approach*3., face, 1.5)+rod(face, face+approach*1.5, 2.75)
                        tool = rod(face+approach*1.5, face+approach*18., 2.)
                        contacts = {name: volume for name, obstacle in blockers.items()
                                    if (volume := float(A.I.intersection_volume(screw+tool, obstacle))) > .001}
                        stock = rod(face-approach*3.5, face, 2.25)
                        coverage = float(A.I.intersection_volume(stock, candidate))/stock.volume()
                        adapter_hits = adapter.ray_cast(face+approach*2., face-approach*20.)
                        adapter_depths = [float((face-np.array(hit.position))@approach)
                                          for hit in adapter_hits]
                        bearing_band = (rod(face-approach*.6, face-approach*.4, 2.75)
                                        -rod(face-approach*.7, face-approach*.3, 1.65))
                        bearing_fraction = float(A.I.intersection_volume(
                            bearing_band, candidate))/bearing_band.volume()
                        blind_trials.append(dict(shift_mm=shift, azimuth_deg=degrees, tilt_deg=tilt,
                                                 face_mm=face.tolist(), axis=approach.tolist(),
                                                 centreline_depth_mm=depth, stock_fraction=coverage,
                                                 adapter_crossings_below_entry_mm=adapter_depths,
                                                 bearing_stock_fraction_at_0_5mm=bearing_fraction,
                                                 hardware_and_tool_hits=contacts,
                                                 tool_cap_overlap_mm3=float(A.I.intersection_volume(tool, candidate))))
            viable = sorted((item for item in blind_trials if not item['hardware_and_tool_hits']
                             and item['tool_cap_overlap_mm3'] <= .001),
                            key=lambda item: -item['stock_fraction'])
            for item in viable:
                face, axis = np.array(item['face_mm']), np.array(item['axis'])
                for radius in (2.5, 3., 3.5):
                    ear = rod(face+axis*.04, face+axis*2.4, radius)
                    hole = rod(face-axis*.2, face+axis*2.7, 1.65)
                    ear = ear-hole
                    merged = (adapter+ear).as_original().simplify(1e-6)
                    added = merged-adapter
                    obstacles = {**blockers, 'trimmed_cap': candidate}
                    collisions = {key: vol for key, part in obstacles.items()
                                  if (vol := float(A.I.intersection_volume(added, part))) > .001}
                    screw_head = rod(face+axis*2.4, face+axis*3.9, 2.5)
                    access = rod(face+axis*3.9, face+axis*20., 2.)
                    access_hits = {key: vol for key, part in {**obstacles, 'adapter': merged}.items()
                                   if (vol := float(A.I.intersection_volume(screw_head+access, part))) > .001}
                    result = dict(shift_mm=item['shift_mm'], azimuth_deg=item['azimuth_deg'],
                                           tilt_deg=item['tilt_deg'], ear_radius_mm=radius,
                                           union_components=sum(part.volume() > 1e-8 for part in merged.decompose()),
                                           existing_adapter_overlap_mm3=float(A.I.intersection_volume(ear, adapter)),
                                           added_stock_hits=collisions, head_and_access_hits=access_hits)
                    ear_trials.append(result)
                    # Recess the head through the adapter and trim the mating
                    # cap only where the new ear occupies it.
                    bore = rod(face-axis*3.2, face+axis*2.41, 1.65)
                    recess = rod(face+axis*2.4, face+axis*22., 2.65)
                    matched_adapter = (merged-hole-recess).as_original().simplify(1e-6)
                    matched_cap = (candidate-ear-recess-bore).as_original().simplify(1e-6)
                    stock = matched_adapter-adapter
                    band = (rod(face+axis*2.2, face+axis*2.4, 2.5)
                            -rod(face+axis*2., face+axis*2.5, 1.65))
                    lost_cap = candidate-matched_cap
                    result['matched_interface'] = dict(
                        bearing_coverage=float(A.I.intersection_volume(band, matched_adapter))/band.volume(),
                        adapter_components=sum(part.volume() > 1e-8 for part in matched_adapter.decompose()),
                        cap_components=sum(part.volume() > 1e-8 for part in matched_cap.decompose()),
                        cap_stock_removed_mm3=float(lost_cap.volume()),
                        adapter_stock_removed_mm3=float((adapter-matched_adapter).volume()),
                        cap_wire_loss_mm3=float(A.I.intersection_volume(lost_cap, wire)),
                        cap_muscle_loss_mm3=float(A.I.intersection_volume(lost_cap, extensor)),
                        residual_added_hits={key: vol for key, part in {**blockers, 'cap': matched_cap}.items()
                                             if (vol := float(A.I.intersection_volume(stock, part))) > .001},
                        access_overlap_mm3=float(A.I.intersection_volume(screw_head+access,
                                                                         matched_adapter+matched_cap)))
                    check = result['matched_interface']
                    if (not exported_interface and check['adapter_components'] == check['cap_components'] == 1
                            and check['bearing_coverage'] >= .99 and not check['residual_added_hits']
                            and max(check['cap_wire_loss_mm3'], check['cap_muscle_loss_mm3'],
                                    check['access_overlap_mm3']) <= .001):
                        try:
                            _, cap_source = export(matched_cap, folder/'anatomical_elbow_matched_cap.obj', face)
                            _, adapter_source = export(
                                matched_adapter, folder/'anatomical_elbow_matched_adapter.obj', face)
                        except (ValueError, AssertionError) as error:
                            check['export_error'] = str(error)
                        else:
                            check['sources'] = dict(cap=cap_source, adapter=adapter_source)
                            check['remaining_adapter_deltoid_overlap_mm3'] = float(
                                A.I.intersection_volume(matched_adapter, deltoid))
                            exported_interface = True
        if (not hits and row['positive_components'] == 1
                and max(row['remaining_overlap_mm3'], row['muscle_interface_loss_mm3'],
                        row['wire_region_loss_mm3']) <= .001):
            origin = np.array(cap.bounding_box()).reshape(2, 3).mean(axis=0)
            _, row['source'] = export(candidate, folder/f'original_extensor_cap_relief_{margin:g}.obj', origin)
    save(folder/'anatomical_elbow_cap_relief.json', dict(
        baseline_sha256=sha(baseline_file), shoulder_clearance_sha256=sha(shoulder_file), trials=rows,
        left_fastener_relocation_trials=relocation,
        angled_blind_fastener_trials=blind_trials,
        adapter_ear_trials=ear_trials,
        scope='Neutral cap relief. Original elbow muscles, attachment points and wire route fixed.'))
    print(rows)


if __name__ == '__main__':
    main()
