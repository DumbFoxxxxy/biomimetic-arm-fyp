"""Recheck inherited hollow guides against the newly unified bones and joint context."""
from itertools import combinations
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component
from build_channel_first_inventory import provenance


def main():
    folder = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    skeleton_file = folder/'unified_skeleton/assembly.json'
    guide_file = folder/'minimal_shoulder_conduits.json'
    skeleton, guides = read(skeleton_file), read(guide_file)
    context_file = folder/'endpoint_neutral_clearance.json'
    context = {**skeleton['sources'], **skeleton['joint_context_sources'],
               **read(context_file)['cervical_context']}
    solids = {name: load_component(source) for name, source in context.items()}
    rows = []
    for guide in guides['conduits']:
        material = load_component(guide['source'])
        wire = A.D.outer_envelope(np.array(guide['path_source_mm']), .3)
        hits, wire_hits = {}, {}
        gap, nearest = 2., None
        for name, solid in solids.items():
            volume = float(A.I.intersection_volume(material, solid))
            if volume > .001:
                hits[name] = volume
            volume = float(A.I.intersection_volume(wire, solid))
            if volume > .001:
                wire_hits[name] = volume
            value = A.I.gap(material, solid, gap)
            if value < gap:
                gap, nearest = value, name
        rows.append(dict(id=guide['id'], material_overlap_mm3=hits, wire_overlap_mm3=wire_hits,
                         material_gap_mm_capped2=gap, nearest_part=nearest))
        print(guide['id'], 'material', hits, 'wire', wire_hits, flush=True)
    rigid_rows = []
    bones = {name: solid for name, solid in solids.items() if name.startswith('bone:')}
    for (left,a),(right,b) in combinations(bones.items(), 2):
        volume = float(A.I.intersection_volume(a,b))
        rigid_rows.append(dict(left=left, right=right, overlap_mm3=volume))
    save(folder/'unified_skeleton/channel_screen.json', dict(
        skeleton=provenance(skeleton_file), guides=provenance(guide_file),
        cervical_context=provenance(context_file), checked_obstacles=len(solids),
        channels=rows, bone_pair_checks=rigid_rows,
        scope='Neutral exported guide material and nominal wire against unified anatomy, cervical context and listed GH/elbow joint parts. End attachments absent; no moving or soft contact proof. Prior dynamic results do not transfer.'))


if __name__ == '__main__':
    main()
