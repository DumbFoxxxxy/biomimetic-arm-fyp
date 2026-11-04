"""Recover collision solids from the two existing elbow TPU meshes."""
from pathlib import Path

import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, mesh_geom, read, save, sha


def main():
    folder = ROOT / '02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus'
    manifest = folder / 'assembly_check.json'
    model = read(manifest)
    mf = A.I._manifold()
    rows = []
    for key in ['tpu:elbow_flexor', 'tpu:elbow_extensor']:
        source = model['sources'][key]
        path = ROOT/source['file']
        assert sha(path) == source['sha256']
        vertices, faces = mesh_geom.load_obj(str(path))
        origin = np.array(source['translation_source_mm'])
        raw = mf.Manifold(mf.Mesh64(np.array(vertices+origin, dtype=np.float64, order='C'),
                                    np.array(faces, dtype=np.uint64, order='C')))
        _, boundary, nonmanifold = mesh_geom.edge_census(faces)
        world_vertices, inverse = np.unique(vertices+origin, axis=0, return_inverse=True)
        indexed = inverse[faces]
        edges = np.sort(np.vstack([indexed[:, [0, 1]], indexed[:, [1, 2]], indexed[:, [2, 0]]]), axis=1)
        unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
        bad_edges = world_vertices[unique_edges[counts > 2]]
        row = dict(key=key, original_source=source, original_boundary_edges=boundary,
                   original_nonmanifold_edges=nonmanifold, manifold_status=str(raw.status()), trials=[],
                   exact_weld_nonmanifold_edges=bad_edges.tolist())
        rows.append(row)
        if raw.is_empty() or str(raw.status()).split('.')[-1] != 'NoError':
            continue
        for tolerance in [0., 1e-6, 1e-5, 1e-4, .001]:
            trial = raw.as_original().simplify(tolerance) if tolerance else raw
            added = max(0., float((trial-raw).volume()))
            removed = max(0., float((raw-trial).volume()))
            record = dict(tolerance_mm=tolerance, added_mm3=added, removed_mm3=removed)
            row['trials'].append(record)
            if added+removed > .01:
                continue
            try:
                _, replacement = export(trial, folder/(key.split(':')[1]+'_collision_candidate.obj'), origin)
            except (AssertionError, ValueError) as error:
                record['error'] = str(error) or type(error).__name__
                continue
            row['candidate_source'] = replacement
            break
        if 'candidate_source' in row or len(bad_edges) == 0:
            continue
        design_file = path.parent/'design.json'
        design = read(design_file)
        route = np.array(design['wire_polyline_display_mm']) + origin
        lumen = A.D.outer_envelope(route, design['tube_inner_radius_mm'])
        row['channel_source_sha256'] = sha(design_file)
        row['edge_patch_trials'] = []
        for radius in [.005, .015, .03]:
            patches = []
            for start, end in bad_edges:
                if np.linalg.norm(end-start) > 1e-8:
                    patches.append(A.D.outer_envelope(np.stack([start, end]), radius))
                else:
                    patches.append(mf.Manifold.sphere(radius, 12).translate(start))
            stock = mf.Manifold.batch_boolean(patches, mf.OpType.Add)
            added = stock-raw
            added_mm3 = max(0., float(added.volume()))
            lumen_overlap = float(A.I.intersection_volume(added, lumen))
            record = dict(radius_mm=radius, added_mm3=added_mm3, lumen_overlap_mm3=lumen_overlap)
            row['edge_patch_trials'].append(record)
            if added_mm3 > .01 or lumen_overlap > .000001:
                continue
            trial = (raw+stock).as_original().simplify(1e-6)
            try:
                _, replacement = export(trial, folder/(key.split(':')[1]+'_collision_candidate.obj'), origin)
            except (AssertionError, ValueError) as error:
                record['error'] = str(error) or type(error).__name__
                continue
            row['candidate_source'] = replacement
            break
    save(folder/'elbow_tpu_readback_repair.json', dict(
        assembly_sha256=sha(manifest), script_sha256=sha(Path(__file__)), rows=rows,
        maximum_material_change_mm3=.01, manufacturing_verified=False))
    print([(row['key'], row['manifold_status'], row.get('candidate_source', {}).get('file'),
            row['trials']) for row in rows], flush=True)


if __name__ == '__main__':
    main()
