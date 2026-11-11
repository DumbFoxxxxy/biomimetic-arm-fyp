"""Locate force-proxy endpoints on allowed actual bones; never freeze drill sites."""
import sys
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save,load_component,mesh_geom
from build_channel_first_inventory import provenance


def main():
    sys.path.insert(0,str(ROOT/'mujoco_arm/anatomical_concept'))
    from extract_attachment_regions import projector
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first';folder=base/'unified_SC_AC'
    endpoint_file=folder/'arm_girdle_equivalent_drive_geometry_mass_check.json';endpoint=read(endpoint_file)
    chain=read(ROOT/endpoint['chain_check']['file']);assembly_file=ROOT/chain['reference']['file']
    assert provenance(assembly_file)==chain['reference']
    assembly=read(assembly_file);sources=dict(assembly['sources'])
    regions=read(base/'thorax_arm_reference.json')['thorax_regions']
    candidates_file=base/'arm_girdle_endpoint_candidates.json'
    candidates={r['id']:r for r in read(candidates_file)['routes']}
    neck_file=ROOT/'02_人体骨骼晶格肌肉方案/输入/骨骼/支撑参考/sources.json'
    neck_models=[r for r in read(neck_file)['models'] if 'cervical_spine' in r['categories']]
    neck={}
    for row in neck_models:
        path=neck_file.parent/f"{row['id']}.obj"
        assert provenance(path)['sha256']==row['sha256']
        key=f"cervical_reference:{row['id']}"
        sources[key]=dict(file=path.relative_to(ROOT).as_posix(),sha256=row['sha256'],translation_source_mm=[0.,0.,0.])
        neck[key]=row['label']
    groups={}
    for body in ('clavicle','scapula','humerus','forearm'):groups[body]=[f'bone:{body}']
    for category in ('ribs','thoracic_spine','sternum'):
        groups[category]=[key for key,r in regions.items() if category in r['categories'] and (category!='ribs' or r['name'].startswith('right'))]
    groups['first_rib']=[key for key in groups['ribs'] if regions[key]['name']=='right first rib']
    groups['cervical_spine']=list(neck)
    definitions=read(ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/joint_model/definition.json')['centres_source_mm']
    cache={}
    def surface(key):
        if key not in cache:
            solid=load_component(sources[key]);v,f=mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64());cache[key]=projector(v,f)
        return cache[key]
    rows=[]
    for route in endpoint['channels']:
        records=[]
        for role,body,point in zip(('origin','distal'),(route['origin_body'],route['distal_body']),route['endpoints_source_mm']):
            structure=candidates[route['id']]['candidate_origin_structure'] if body=='thorax' else body
            assert groups[structure],structure
            options=[]
            for key in groups[structure]:
                projected,distance,_=surface(key)(np.array([point]));options.append((float(distance[0]),key,projected[0]))
            distance,key,q=min(options,key=lambda r:r[0])
            zones=[]
            for joint,radius in [('SC',36.),('AC',36.),('GH',26.),('elbow',45.)]:
                if np.linalg.norm(q-np.array(definitions[joint]))<radius:zones.append(joint)
            records.append(dict(role=role,body=body,allowed_structure=structure,
                proxy_point_source_mm=point,nearest_bone_key=key,
                bone_name=neck.get(key,regions.get(key,{}).get('name',key)),
                surface_point_source_mm=q.tolist(),gap_mm=distance,
                bone_present_in_selected_assembly=key in assembly['sources'],
                engineering_joint_regions=zones,
                status='missing_physical_support' if key not in assembly['sources'] else
                       'joint_interface_anchor_review_required' if zones else
                       'adapter_required' if distance>.5 else 'surface_candidate_hardware_pending',
                frozen_anchor=False))
        rows.append(dict(id=route['id'],endpoints=records,drive_end=None,supply_channel=None,
                         reaction_seat=None,terminal_fixation=None,
                         note='Actual motor, sheath support and terminal hardware remain to be designed; projection is not a selected anchor.'))
    save(folder/'anchor_support_audit.json',dict(endpoint_model=provenance(endpoint_file),
        assembly=provenance(assembly_file),anatomical_candidates=provenance(candidates_file),
        cervical_inventory=provenance(neck_file),reference_cervical_sources={k:sources[k] for k in neck},
        generator=provenance(ROOT/'tools/audit_channel_first_anchor_supports.py'),routes=rows,
        scope='54 force-proxy endpoints projected only onto anatomically allowed bone groups. Nearest projection does not choose a drill hole, terminal, safe bone wall or actual origin/insertion. Missing cervical supports and joint-region conflicts remain explicit. No endpoints or source geometry changed.'))
    missing=[(r['id'],p['bone_name']) for r in rows for p in r['endpoints'] if not p['bone_present_in_selected_assembly']]
    print('Checked',len(rows)*2,'endpoints; missing physical bone supports',missing,flush=True)
    print('Joint-region reviews',sum(bool(p['engineering_joint_regions']) for r in rows for p in r['endpoints']),flush=True)


if __name__=='__main__':main()
