"""Apply compact AC engineering deltas while retaining current SC and GH stock."""
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import A, ROOT, read, save, load_component, export
from build_channel_first_inventory import provenance


def main():
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=folder/'unified_SC/assembly.json';parent=read(parent_file)
    legacy=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/SC_connection'
    design_file=legacy/'AC_compact/assembly.json';design=read(design_file)
    centre=np.array(design['centre_source_mm']);dest=folder/'unified_SC_AC';dest.mkdir(exist_ok=True)
    sources,joints=dict(parent['sources']),dict(parent['joint_context_sources'])
    def imported(file,translation):
        meta=dict(file=file.relative_to(ROOT).as_posix(),sha256=provenance(file)['sha256'],translation_source_mm=translation)
        return load_component(meta),meta
    edits=[]
    for key,name,original_file,shift in [
        ('bone:clavicle','clavicle_SC_AC_balls',legacy/'clavicle_SC_ball.obj',[-140.,-70.,1100.]),
        ('bone:scapula','scapula_AC_socket',ROOT/'02_人体骨骼晶格肌肉方案/输入/骨骼/肩胛骨.obj',[0.,0.,0.])]:
        old,old_meta=imported(original_file,shift)
        new,new_meta=imported(legacy/'AC_compact'/f'{name}.obj',[-140.,-70.,1100.])
        region=A.I._manifold().Manifold.sphere(45,64).translate(centre)
        added=(new-old)^region;removed=(old-new)^region
        current=load_component(sources[key])
        candidate=(current-removed+added).as_original().simplify(1e-5)
        if key=='bone:scapula':
            definition=read(ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/definition.json')
            gh=np.array(definition['centres_source_mm']['GH'])
            # Transfer a volumetric patch rather than sub-tolerance differences
            # between two independently simplified copies of the bone skin.
            protect=A.I._manifold().Manifold.sphere(26.,64).translate(gh)
            inner=A.I._manifold().Manifold.sphere(35.,64).translate(centre)-A.I._manifold().Manifold.sphere(27.,64).translate(gh)
            outer=A.I._manifold().Manifold.sphere(36.,64).translate(centre)-protect
            candidate=((current-inner)+(new^outer)).as_original().simplify(1e-5)
        attempts=[]
        for tolerance in (1e-4,.001,.01):
            prepared=candidate.as_original().simplify(tolerance)
            try:
                actual,meta=export(prepared,dest/f'{name}.obj',centre)
            except (AssertionError,ValueError) as error:
                attempts.append(dict(tolerance_mm=tolerance,error=str(error) or 'exact mesh topology rejected'))
                continue
            change=max(0.,float((actual-candidate).volume()))+max(0.,float((candidate-actual).volume()))
            attempts.append(dict(tolerance_mm=tolerance,symmetric_material_change_mm3=change))
            if change<.1:break
        else:
            (dest/f'{name}.obj').unlink(missing_ok=True)
            try:
                actual,meta=export(candidate.as_original().simplify(.001),dest/f'{name}.stl',centre)
                change=max(0.,float((actual-candidate).volume()))+max(0.,float((candidate-actual).volume()))
                attempts.append(dict(format='STL',symmetric_material_change_mm3=change,
                                     import_note='float32 export and checked standard mesh import'))
                assert change<.1
            except (AssertionError,ValueError) as error:
                attempts.append(dict(format='STL',error=str(error) or 'standard mesh import rejected'))
                save(dest/'generation_failure.json',dict(key=key,attempts=attempts,status='no_approved_export'))
                (dest/f'{name}.stl').unlink(missing_ok=True)
                raise ValueError('AC delta export rejected; see generation_failure.json') from error
        sources[key]=meta
        protection={}
        if key=='bone:scapula':
            protection=dict(GH_protection_radius_mm=26.,
                GH_protected_material_change_mm3=max(0.,float(((actual-current)^protect).volume()))+max(0.,float(((current-actual)^protect).volume())),
                method='35/36 mm AC volumetric patch with 27/26 mm GH protection and overlap band')
        edits.append(dict(key=key,old_design_source=old_meta,new_design_source=new_meta,
            transferred_addition_mm3=float(added.volume()),transferred_removal_mm3=float(removed.volume()),
            actual_addition_mm3=float((actual-current).volume()),actual_removal_mm3=float((current-actual).volume()),export_attempts=attempts,protection=protection))
    for name in design['components']:
        if name in ('clavicle_SC_AC_balls','scapula_AC_socket'):continue
        solid,ref=imported(legacy/'AC_compact'/f'{name}.obj',[-140.,-70.,1100.])
        _,meta=export(solid,dest/f'{name}.obj',centre);joints[f'AC:{name}']=meta
    solids={k:load_component(v) for k,v in {**sources,**joints}.items()}
    def owner(k):
        if k.startswith('bone:'):return k[5:]
        if k.startswith('SC:'):return 'thorax'
        if k=='elbow:elbow_bearing_bush':return 'forearm'
        if k.startswith('elbow:'):return 'humerus'
        if k in joints:return 'scapula'
        return 'thorax'
    hits=[];count=0
    for a,b in combinations(solids,2):
        if owner(a)==owner(b) and not(a.startswith('AC:') or b.startswith('AC:')):continue
        count+=1;volume=float(A.I.intersection_volume(solids[a],solids[b]))
        if volume>.001:hits.append(dict(left=a,right=b,overlap_mm3=volume))
    save(dest/'assembly.json',dict(sources=sources,joint_context_sources=joints,parent=provenance(parent_file),
        AC_design=provenance(design_file),local_changes=edits,AC_centre_source_mm=centre.tolist(),
        checked_pairs=count,baseline_interferences=hits,nominal_ball_radius_mm=6.,nominal_socket_radius_mm=6.35,
        status='AC_delta_candidate_pending_fit_and_motion_approval',
        scope='Local compact AC geometry delta, including connected supports, split caps and fastener envelopes. No bone or joint centre moved. Interferences retained, not filtered out. No assembly/strength/driven-range approval.'))
    save(dest/'joint_context_scene.json',dict(sources={**sources,**joints},status='static_AC_delta_candidate'))
    print('AC candidate',len(solids),'objects;',count,'pairs; hits',hits,flush=True)


if __name__=='__main__':main()
