"""Evaluate nonnegative scapulothoracic support forces at geometry-derived gaps."""
import sys
import argparse

import numpy as np
from scipy.optimize import linprog, lsq_linear

from build_s1_teres_origin_mount import ROOT, mesh_geom, read, save
from build_channel_first_inventory import provenance


def main():
    import mujoco

    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT/'mujoco_arm/anatomical_concept'))
    from project_runtime import load_mujoco_model
    from extract_attachment_regions import projector

    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--unified-sc-ac',action='store_true')
    args=parser.parse_args()
    root_folder=folder
    if args.unified_sc_ac:folder=folder/'unified_SC_AC'
    model_report_file = folder/'arm_girdle_equivalent_drive_geometry_mass_check.json'
    model_report = read(model_report_file)
    model_file = ROOT/model_report['generated_model']['file']
    assert provenance(model_file) == model_report['generated_model']
    reference_file = folder/'thorax_arm_reference.json'
    if args.unified_sc_ac:
        chain=read(ROOT/model_report['chain_check']['file'])
        reference_file=ROOT/chain['reference']['file']
        assert provenance(reference_file)==chain['reference']
    reference = read(reference_file)
    thorax_regions=read(root_folder/'thorax_arm_reference.json')['thorax_regions']
    scapula_source = reference['sources']['bone:scapula']
    path = ROOT/scapula_source['file']
    assert provenance(path)['sha256'] == scapula_source['sha256']
    vertices, faces = (mesh_geom.load_obj(str(path)) if path.suffix == '.obj'
                       else mesh_geom.read_stl(str(path)))
    vertices += np.array(scapula_source.get('translation_source_mm') or [0., 0., 0.])
    project_scapula = projector(vertices, faces)
    rib_vertices, rib_faces = [], []
    count = 0
    for key, region in thorax_regions.items():
        if 'ribs' not in region['categories'] or not region['name'].startswith('right'):
            continue
        source = reference['sources'][key]
        file = ROOT/source['file']
        assert provenance(file)['sha256'] == source['sha256']
        v, f = mesh_geom.load_obj(str(file))
        rib_vertices.append(v)
        rib_faces.append(f+count)
        count += len(v)
    project_ribs = projector(np.vstack(rib_vertices), np.vstack(rib_faces))
    patch_file = ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/attachment_regions/FJ1459_insertion_scapula.npz'
    with np.load(patch_file) as patch:
        points = patch['bone_points_source_mm']
    points = points[np.linspace(0, len(points)-1, min(240, len(points))).astype(int)]
    on_scapula, _, _ = project_scapula(points)
    on_ribs, gaps, _ = project_ribs(on_scapula)
    valid = np.flatnonzero((gaps >= .5) & (gaps <= 15.))
    # Sample broad attachment extent before inspecting required force.
    indices = []
    while len(indices) < 12 and len(valid):
        if not indices:
            chosen = valid[np.argmin(gaps[valid])]
        else:
            distances = np.linalg.norm(on_scapula[valid, None]-on_scapula[indices], axis=2).min(axis=1)
            if distances.max() < 12.:
                break
            chosen = valid[np.argmax(distances)]
        indices.append(int(chosen))
        valid = valid[valid != chosen]
    model = load_mujoco_model(model_file)
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    matrix = np.zeros((model.nv, len(model_report['channels'])))
    for block in model_report['neutral_joint_blocks'].values():
        matrix[block['velocity_dof_indices']] = np.array(block['moment_arms_mm']).T*.001
    support_columns = []
    candidates = []
    body_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, 'scapula')
    for index in indices:
        point, rib, gap = on_scapula[index], on_ribs[index], gaps[index]
        normal = (point-rib)/gap
        jacobian = np.zeros((3, model.nv))
        rotation = np.zeros_like(jacobian)
        mujoco.mj_jac(model, data, jacobian, rotation, point*.001, body_id)
        support_columns.append(jacobian.T@normal)
        candidates.append(dict(scapular_surface_source_mm=point.tolist(), rib_surface_source_mm=rib.tolist(),
                               normal_force_on_scapula_source=normal.tolist(), gap_mm=float(gap),
                               contact='candidate pad fills this gap; physical attachment and clearance pending'))
    target = data.qfrc_bias-data.qfrc_passive
    combined = np.c_[matrix, np.array(support_columns).T]
    fits = []
    for support_limit in (20., 60., 120.):
        upper = np.r_[np.full(matrix.shape[1], 60.), np.full(len(candidates), support_limit)]
        fit = lsq_linear(combined, target, bounds=(np.zeros(len(upper)), upper), tol=1e-10, max_iter=1000)
        force = fit.x[matrix.shape[1]:]
        fits.append(dict(support_force_limit_N=support_limit,
                         cable_tensions_N=fit.x[:matrix.shape[1]].tolist(), support_forces_N=force.tolist(),
                         residual_norm_Nm=float(np.linalg.norm(fit.fun)),
                         residual_by_joint_Nm={name: fit.fun[block['velocity_dof_indices']].tolist()
                                               for name, block in model_report['neutral_joint_blocks'].items()}))
    first = fits[0]
    initial_forces = np.r_[first['cable_tensions_N'], first['support_forces_N']]
    sparse = linprog(np.ones(combined.shape[1]), A_eq=combined, b_eq=combined@initial_forces,
                     bounds=[(0., 60.)]*matrix.shape[1]+[(0., 20.)]*len(candidates), method='highs')
    compact = None
    if sparse.success:
        required = sparse.x[matrix.shape[1]:]
        compact = dict(
            candidate_indices=np.flatnonzero(required > 1e-5).tolist(),
            cable_tensions_N=sparse.x[:matrix.shape[1]].tolist(),
            support_forces_N=required.tolist(),
            residual_norm_Nm=float(np.linalg.norm(combined@sparse.x-target)),
            objective='minimum total force at the same achievable reaction vector',
            required_pad_compression_mm_at_2N_per_mm=(required/2.).tolist(),
            pad_stiffness_assumption_N_per_mm=2.)
    save(folder/'scapulothoracic_support_direction_screen.json', dict(
        endpoint_model=provenance(model_report_file), anatomy=provenance(reference_file),
        anatomical_patch=provenance(patch_file), candidates=candidates, fits=fits,
        compact_force_solution=compact,
        assumptions=dict(min_gap_mm=.5, max_gap_mm=15., max_candidates=12,
                         nonnegative_normal_forces_only=True, tangential_force=False),
        scope='Neutral reaction-direction feasibility only. Required pads, mounts, compression law, '
              'contact activation and moving contact are unbuilt. No physical equilibrium or motion pass.'))
    print('candidate normal supports', len(candidates),
          [(row['support_force_limit_N'], row['residual_norm_Nm']) for row in fits])


if __name__ == '__main__':
    main()
