"""Restore four explicit M3 bosses/bores for unchanged original elbow supports."""
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,export,rod,mesh_geom
from build_channel_first_inventory import provenance


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'elbow_interface_integration/assembly.json';parent=read(parent_file)
    design_file=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/elbow_muscles/mounted_elbow_candidate/assembly.json'
    design=read(design_file);sources=dict(parent['sources']);bone=load_component(sources['bone:humerus'])
    bosses=[];bores=[];rows=[]
    for muscle in ('FJ1486','FJ1477'):
        placement=next(r for r in design['placements'] if r['id']==muscle)
        for i,anchor in enumerate(placement['bone_anchors']):
            key=f'elbow:{muscle}_bone_{i}_bolt';source=sources[key]
            translation=np.array(source['translation_source_mm']);normal=np.array(anchor['bolt_axis'])
            front=np.array(anchor['front_source_bone_mm'])+translation
            back=np.array(anchor['back_source_bone_mm'])+translation
            qfront=np.array(anchor['front_boss_face_mm'])+translation
            qback=np.array(anchor['back_boss_face_mm'])+translation
            bolt=load_component(source);v,_=mesh_geom.manifold_mesh_to_arrays(bolt.to_mesh64())
            offset=v-front;radial=np.linalg.norm(offset-(offset@normal)[:,None]*normal,axis=1)
            assert radial.max()<2.751,'Mount metadata axis does not match actual bolt'
            bosses += [rod(front-normal*3.,qfront,4.5),rod(back+normal*3.,qback,4.5)]
            bores.append(rod(qback-normal*20.,qfront+normal*20.,1.7))
            rows.append(dict(muscle=muscle,anchor=i,bolt_source=source,axis_point_source_mm=front.tolist(),
                axis_source=normal.tolist(),front_seat_source_mm=qfront.tolist(),back_seat_source_mm=qback.tolist(),
                boss_radius_mm=4.5,bore_radius_mm=1.7))
    mf=A.I._manifold();candidate=(mf.Manifold.batch_boolean([bone,*bosses],mf.OpType.Add)-mf.Manifold.batch_boolean(bores,mf.OpType.Add)).as_original().simplify(1e-6)
    dest=base/'elbow_mount_candidate';dest.mkdir(exist_ok=True)
    centre=np.mean([r['axis_point_source_mm'] for r in rows],axis=0)
    actual,source=export(candidate,dest/'humerus_with_mount_seats.obj',centre)
    sources['bone:humerus']=source
    definition=read(ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/definition.json')['centres_source_mm']
    protection={}
    for joint in ('GH','elbow'):
        region=mf.Manifold.sphere(45.,64).translate(definition[joint])
        change=max(0.,float(((actual-bone)^region).volume()))+max(0.,float(((bone-actual)^region).volume()))
        assert change<.001;protection[joint]=change
    hits=[]
    for key,meta in sources.items():
        if key=='bone:humerus':continue
        volume=float(A.I.intersection_volume(actual,load_component(meta)))
        if volume>.001:hits.append(dict(target=key,overlap_mm3=volume))
    save(dest/'assembly.json',dict(sources=sources,display_only_TPU_sources=parent['display_only_TPU_sources'],
        legacy_owners=parent['legacy_owners'],muscles=parent['muscles'],parent=provenance(parent_file),
        mounting_design=provenance(design_file),generator=provenance(ROOT/'tools/build_channel_first_elbow_mount_seats.py'),
        mount_seats=rows,added_mm3=float((actual-bone).volume()),removed_mm3=float((bone-actual).volume()),
        protected_joint_material_change_mm3=protection,humerus_interface_hits=hits,
        status='restored_mount_candidate_pending_motion_and_assembly',
        scope='Four original M3 axes, 4.5mm bearing bosses, 1.7mm clearance bores. Original elbow caps, lattice, cable routes and attachment locations unchanged. Threads, preload, strength, tools and passive swivel mechanics remain unverified.'))
    save(dest/'view_scene.json',dict(sources={**sources,**parent['display_only_TPU_sources']},status='restored_original_elbow_mount_candidate'))
    print('Four mounts restored; hits',hits,'protected GH/elbow changes',protection,flush=True)


if __name__=='__main__':main()
