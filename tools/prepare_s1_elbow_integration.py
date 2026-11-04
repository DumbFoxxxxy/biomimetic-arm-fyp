"""Prepare a LOCAL spatial overlay and shared-humerus integration audit.

Not an assembled/frozen design and never a publishing command. Keep the live
GitHub viewer unchanged until S1 and its common-bone interfaces are finalized.
"""
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from project_runtime import project_path
from build_s1_teres_origin_mount import ROOT,A,read,save,sha,load_component

DEST=ROOT/'02_人体骨骼晶格肌肉方案/输出/integration_review/s1_elbow_preview'
ELBOW=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers/continuous_channel_candidate/neutral_view.xml'


def main():
    DEST.mkdir(parents=True,exist_ok=True)
    shoulder_file=A.DESIGN/'s1_posterior_split_clamp/assembly_check.json'
    shoulder=read(shoulder_file);proof_file=shoulder_file.parent/'verification_checks.json';proof=read(proof_file)
    assert proof['assembly_sha256']==sha(shoulder_file) and proof['local_rigid_approaches_and_seats_pass']
    for row in shoulder['sources'].values():assert sha(ROOT/row['file'])==row['sha256']
    tree=ET.parse(ELBOW).getroot()
    assert not tree.findall('.//body'),'Nested elbow scene needs explicit transform handling'
    meshes={row.get('name'):row for row in tree.findall('./asset/mesh')}
    rows={};old_shoulder=[];elbow_humerus=None
    for geom in tree.findall('./worldbody/geom'):
        name=geom.get('name');asset=meshes[geom.get('mesh')]
        assert geom.get('type')=='mesh' and not any(geom.get(k) for k in ['quat','euler','axisangle','xyaxes','zaxis'])
        scale=np.fromstring(asset.get('scale','1 1 1'),sep=' ')
        assert np.allclose(scale,[.001,.001,.001],rtol=0,atol=1e-12),'Unexpected elbow mesh scale'
        path=Path(project_path(asset.get('file')))
        assert path.is_file() and path.is_relative_to(ROOT)
        offset=A.frames.SHIFT+np.fromstring(geom.get('pos','0 0 0'),sep=' ')*1000.
        row=dict(file=path.relative_to(ROOT).as_posix(),sha256=sha(path),translation_source_mm=offset.tolist(),
                 source_elbow_geom=name,mesh_scale_m_per_file_unit=scale.tolist())
        if name=='FJ3368':elbow_humerus=row;old_shoulder.append(name);continue
        if name in ['FJ3362','FJ3384','FJ3290'] or name.startswith(('SC_','AC_','GH_')):
            old_shoulder.append(name);continue
        if name=='FJ3391':key='bone:forearm'
        elif name=='FJ1486_TPU_muscle':key='tpu:elbow_flexor'
        elif name=='FJ1477_TPU_muscle':key='tpu:elbow_extensor'
        else:key='elbow:'+name
        assert key not in shoulder['sources'] and key not in rows
        rows[key]=row
    assert elbow_humerus is not None
    sources={**shoulder['sources'],**rows}
    overlay=dict(scope='LOCAL spatial overlay only: S1 humerus retained; elbow humerus NOT merged. Interfaces and inter-module collisions pending. Not a frozen assembly.',
        sources=sources,n_instances=len(sources),source_shoulder=shoulder_file.relative_to(ROOT).as_posix(),
        source_shoulder_sha256=sha(shoulder_file),source_elbow_scene=ELBOW.relative_to(ROOT).as_posix(),source_elbow_scene_sha256=sha(ELBOW),
        shoulder_objects=len(shoulder['sources']),elbow_objects_added=len(rows),elbow_legacy_shoulder_context_omitted=old_shoulder,
        shared_humerus_policy='Only S1 humerus displayed. Elbow humerus kept separately as comparison input; no Boolean union or feature deletion.',
        display_to_source_translation_mm=A.frames.SHIFT.tolist(),
        exact_mesh_coordinates_preserved=True,shared_humerus_integrated=False,combined_assembly_verified=False,
        public_release_ready=False,stage_pass=False)
    save(DEST/'spatial_overlay.json',overlay)
    S=load_component(shoulder['sources']['bone:humerus']);E=load_component(elbow_humerus)
    differences=[]
    for key,solid in [('elbow_only_material',E-S),('S1_only_material',S-E)]:
        differences.append(dict(region=key,volume_mm3=float(solid.volume()),bounds_source_mm=list(solid.bounding_box()),
            components_over_1mm3=[dict(volume_mm3=float(part.volume()),bounds_source_mm=list(part.bounding_box()))
                                   for part in solid.decompose() if abs(part.volume())>1.]))
    fit=[];errors=[]
    for key,row in rows.items():
        try:part=load_component(row)
        except Exception as error:
            errors.append(dict(key=key,error=str(error)));continue
        volume=float(A.I.intersection_volume(S,part));old=float(A.I.intersection_volume(E,part))
        if volume>.001 or old>.001:
            fit.append(dict(key=key,S1_humerus_overlap_mm3=volume,original_elbow_humerus_overlap_mm3=old,
                            increased_over_0p001_mm3=volume>max(.001,old+.001)))
    audit=dict(scope=__doc__,script_sha256=sha(Path(__file__)),overlay_sha256=sha(DEST/'spatial_overlay.json'),
        shoulder_sha256=sha(shoulder_file),shoulder_local_proof_sha256=sha(proof_file),elbow_xml_sha256=sha(ELBOW),
        original_elbow_humerus=elbow_humerus,displayed_S1_humerus=shoulder['sources']['bone:humerus'],
        same_coordinate_frame=True,scale_adjustment_applied=False,shared_humerus_differences=differences,
        elbow_parts_checked_against_shared_humerus=len(rows)-len(errors),readback_errors=errors,humerus_interface_overlaps=fit,
        whole_combined_collision_audit_completed=False,common_humerus_geometry_unified=False,
        S1_freeze_pending=['anterior deltoid bone attachment','middle deltoid bone attachment','shared distal deltoid attachment',
                           'whole shoulder assembly consistency and joint/muscle clearance review'],
        publication_condition='Finalize S1 mechanical geometry; unify the shared humerus and check combined interfaces; label static/predicted/unverified capabilities accurately. No physical printing or full FE prerequisite added.',
        public_release_ready=False,live_pages_assets_modified=False)
    save(DEST/'integration_audit.json',audit)
    print(dict(shoulder=len(shoulder['sources']),elbow_added=len(rows),overlay=len(sources),
        original_elbow_geoms=len(tree.findall('./worldbody/geom')),legacy_context_omitted=len(old_shoulder),
        increased_humerus_overlaps=sum(row['increased_over_0p001_mm3'] for row in fit),readback_errors=len(errors),public_release_ready=False),flush=True)


if __name__=='__main__':main()
