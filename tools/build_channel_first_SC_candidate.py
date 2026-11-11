"""Transfer a real SC ball/socket locally into the unified channel-first bones."""
from itertools import combinations
import numpy as np

from build_s1_teres_origin_mount import A, ROOT, read, save, load_component, export
from build_channel_first_inventory import provenance


def main():
    base = ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file = base/'unified_skeleton/assembly.json'
    parent = read(parent_file)
    legacy = ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/SC_connection'
    legacy_file = legacy/'assembly.json'
    design = read(legacy_file)
    centre = np.array(design['joint_centre_source_mm'])
    dest = base/'unified_SC'
    dest.mkdir(exist_ok=True)
    sources, joints = dict(parent['sources']), dict(parent['joint_context_sources'])
    imported = {}
    for name in design['components']:
        file = legacy/f'{name}.obj'
        source = dict(file=file.relative_to(ROOT).as_posix(),sha256=provenance(file)['sha256'],
                      translation_source_mm=[-140.,-70.,1100.])
        imported[name] = (load_component(source), source)
    deltas = []
    inner = A.I._manifold().Manifold.sphere(35.,64).translate(centre)
    outer = A.I._manifold().Manifold.sphere(36.,64).translate(centre)
    for key,name in [('bone:clavicle','clavicle_SC_ball'),('thorax_reference:FJ3290','manubrium_socket_base')]:
        original = load_component(sources[key])
        candidate = ((original-inner)+(imported[name][0]^outer)).as_original().simplify(1e-6)
        actual, source = export(candidate,dest/f'{name}.obj',centre)
        outside_added = max(0.,float(((actual-original)-outer).volume()))
        outside_removed = max(0.,float(((original-actual)-outer).volume()))
        assert outside_added+outside_removed < .001
        deltas.append(dict(key=key,added_mm3=float((actual-original).volume()),
                           removed_mm3=float((original-actual).volume()),
                           outside_36mm_change_mm3=outside_added+outside_removed,
                           local_design_source=imported[name][1]))
        sources[key] = source
    for name,(solid,reference) in imported.items():
        if name in ('clavicle_SC_ball','manubrium_socket_base'):continue
        _,source=export(solid,dest/f'{name}.obj',centre)
        joints[f'SC:{name}']=source
    all_sources={**sources,**joints}
    solids={name:load_component(source) for name,source in all_sources.items()}
    def owner(key):
        if key.startswith('bone:'):return key[5:]
        if key.startswith('SC:'):return 'thorax'
        if key=='elbow:elbow_bearing_bush':return 'forearm'
        if key.startswith('elbow:'):return 'humerus'
        if key in joints:return 'scapula'
        return 'thorax'
    hits=[];count=0
    for a,b in combinations(solids,2):
        if owner(a)==owner(b) and not (a.startswith('SC:') or b.startswith('SC:')):continue
        count+=1
        volume=float(A.I.intersection_volume(solids[a],solids[b]))
        if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
    report=dict(sources=sources,joint_context_sources=joints,parent=provenance(parent_file),
        SC_design=provenance(legacy_file),SC_centre_source_mm=centre.tolist(),local_changes=deltas,
        nominal_ball_radius_mm=10.,nominal_socket_radius_mm=10.35,
        split_caps=2,nominal_M3_fastener_sets=2,checked_pairs=count,baseline_interferences=hits,
        status='SC_local_candidate_pending_motion_and_mount_validation',
        scope='Locally transferred ball, connected neck, socket support, split retention caps and fastener envelopes. Outside 36 mm sphere material preserved. No joint centre moved, no adjacency pair excluded. Threads, strength, tolerance and physical SC/AC assembly remain unverified.')
    save(dest/'assembly.json',report)
    save(dest/'joint_context_scene.json',dict(sources=all_sources,status=report['status']))
    print('SC candidate',len(all_sources),'objects;',count,'pairs; hits',hits,flush=True)


if __name__=='__main__':main()
