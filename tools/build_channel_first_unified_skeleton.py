"""Combine selected S1 proximal stock with original-elbow distal stock at a verified seam."""
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, export, load_component, read, save
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    dest = folder/'unified_skeleton'
    dest.mkdir(exist_ok=True)
    ref_file = folder/'thorax_arm_reference.json'
    selection_file = ROOT/'design/S1_CURRENT_CANDIDATE.json'
    selection = read(selection_file)
    shoulder_file = ROOT/selection['assembly_report']
    assert provenance(shoulder_file)['sha256'] == selection['assembly_sha256']
    shoulder, reference = read(shoulder_file), read(ref_file)
    elbow_file = ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_shared_humerus/assembly_check.json'
    elbow = read(elbow_file)
    frame_file = A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame = read(frame_file)['frame_source']
    basis, origin = np.array(frame['basis_columns']), np.array(frame['origin_mm'])
    local, world = np.c_[basis.T, -basis.T@origin], np.c_[basis, origin]
    upper_source = dict(shoulder['sources']['bone:humerus'])
    upper_source['translation_source_mm'] = upper_source.get('translation_source_mm') or [0., 0., 0.]
    upper = load_component(upper_source).transform(local)
    lower = load_component(reference['sources']['bone:humerus']).transform(local)
    cube = A.I._manifold().Manifold.cube
    band = cube([400,400,4]).translate([-200,-200,108])
    ua, la = upper^band, lower^band
    seam_difference = float((ua-la).volume()+(la-ua).volume())
    seam_common = float((ua^la).volume())
    assert abs(seam_difference) < .001 and seam_common > 100.
    proximal = cube([400,400,800]).translate([-200,-200,-689])  # ends at station111
    distal = cube([400,400,800]).translate([-200,-200,109])
    candidate = ((upper^proximal)+(lower^distal)).transform(world).as_original().simplify(1e-6)
    actual, source = export(candidate, dest/'humerus.obj', origin)
    posed = actual.transform(local)
    upper_check = cube([400,400,800]).translate([-200,-200,-692])  # ends108
    lower_check = cube([400,400,800]).translate([-200,-200,112])
    def difference(a, b):
        return max(0., float((a-b).volume()))+max(0., float((b-a).volume()))
    preserved_upper = difference(posed^upper_check, upper^upper_check)
    preserved_lower = difference(posed^lower_check, lower^lower_check)
    assert preserved_upper < .001 and preserved_lower < .001
    sources = dict(reference['sources'])
    for name in ('bone:clavicle','bone:scapula'):
        sources[name] = dict(shoulder['sources'][name])
        sources[name]['translation_source_mm'] = sources[name].get('translation_source_mm') or [0.,0.,0.]
    sources['bone:humerus'] = source
    joints = {k:v for k,v in shoulder['sources'].items() if k.startswith(('GH_cap_', 'GH_bolt_', 'GH_nut_', 'GH_top_washer_', 'GH_bottom_washer_', 'GH_keeper:'))}
    joints.update({k:v for k,v in elbow['sources'].items() if k.startswith('elbow:elbow_')})
    fit = []
    for name, metadata in {**{k:v for k,v in sources.items() if k.startswith('bone:') and k!='bone:humerus'}, **joints}.items():
        metadata = dict(metadata)
        metadata['translation_source_mm'] = metadata.get('translation_source_mm') or [0.,0.,0.]
        volume = float(A.I.intersection_volume(actual, load_component(metadata)))
        fit.append(dict(target=name, overlap_mm3=volume))
    save(dest/'assembly.json', dict(sources=sources, joint_context_sources=joints,
        reference=provenance(ref_file), shoulder_selection=provenance(selection_file),
        shoulder_assembly=provenance(shoulder_file), elbow_assembly=provenance(elbow_file),
        seam_frame=provenance(frame_file), generator=provenance(ROOT/'tools/build_channel_first_unified_skeleton.py'),
        seam=dict(station_mm=110., overlap_interval_mm=[109.,111.], check_band_mm=[108.,112.],
                  symmetric_difference_mm3=seam_difference, common_stock_volume_mm3=seam_common),
        upper_preservation_difference_mm3=preserved_upper, lower_preservation_difference_mm3=preserved_lower,
        humerus_interface_checks=fit,
        static_humerus_interface_clear=all(row['overlap_mm3']<=.001 for row in fit),
        status='unified_geometry_candidate_not_motion_or_assembly_approved',
        pending=['SC/AC physical support on thorax and tabletop base', 'all rigid-pair continuous motion',
                 'end fixtures on revised bones', 'all channel results must be recomputed on this reference',
                 'editable engineering CAD and physical assembly'],
        scope='One continuous humerus; selected S1 upper joint and original elbow distal features preserved. No bone or joint transform changed. No backboard, bench clamps or old muscle fixtures counted as installed. Joint context is separately listed; only humerus interfaces checked here.'))
    save(dest/'joint_context_scene.json', dict(sources={**sources,**joints},
        status='static_bones_and_joint_context_not_full_assembly'))
    print('Unified humerus', 'seam difference', seam_difference,
          'upper/lower preservation', preserved_upper, preserved_lower,
          'interface hits', [r for r in fit if r['overlap_mm3']>.001], flush=True)


if __name__ == '__main__':
    main()
