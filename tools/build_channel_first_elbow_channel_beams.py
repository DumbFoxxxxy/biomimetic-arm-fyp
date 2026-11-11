"""Embed compliant supported channel chains into the actual passive-cap model."""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector
from simulate_channel_first_elastic_conduits import build_model


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    parent_file=base/'elbow_passive_model/model_check.json';parent=read(parent_file)
    source_file=ROOT/parent['model']['file'];assert provenance(source_file)==parent['model']
    assembly_file=ROOT/parent['assembly']['file'];assembly=read(assembly_file)
    assert provenance(assembly_file)==parent['assembly']
    root=ET.parse(source_file).getroot();root.find('option').set('timestep','.00002')
    root.find('option').set('jacobian','sparse')
    bodies={b.get('name'):b for b in root.findall('.//body')}
    equality=root.find('equality')
    if equality is None:equality=ET.SubElement(root,'equality')
    rows=[]
    for label,muscle in assembly['muscles'].items():
        head=muscle['anatomical_id'];path=np.array(muscle['wire_path_source_mm'][1:-1])
        beam,spec=build_model(dict(path_source_mm=path.tolist()),15.,5.,.00002)
        subtree=beam.find("./worldbody/body[@name='proximal_clamp']")
        # Prefix all body/joint/site names, preserving rest directions and elasticity.
        for element in subtree.iter():
            if element.get('name'):element.set('name',head+'__'+element.get('name'))
        origin=np.array(spec['source_origin_m'])
        subtree.set('pos',vector(origin-np.array(parent['neutral_origins_m']['humerus'])))
        bodies['humerus'].append(subtree)
        count=spec['segments']
        ET.SubElement(equality,'connect',name=head+'_channel_to_cap',
            site1=f'{head}__node_{count}',site2=head+'_distal',solref='.001 1',solimp='.999 .9999 .0001')
        tendon=root.find(f"./tendon/spatial[@name='{head}']")
        for child in list(tendon):tendon.remove(child)
        tendon.set('stiffness','0');tendon.set('damping','0')
        for i in range(count+1):ET.SubElement(tendon,'site',site=f'{head}__node_{i}')
        rows.append(dict(head=head,path_source_mm=path.tolist(),specification=spec,
            proximal_support='clamped to existing humerus-owned proximal rigid cap frame',
            distal_support='translation connected to original bore on passive distal cap; rotation free',
            intermediate_support='moving compliant beam elements, no world-fixed via points',
            parameter_note='Annular radius1.6/0.8mm reference, E15MPa, not full original lattice stiffness; existing tendon200N/m spring removed for these channels to avoid double-counting'))
    folder=base/'elbow_channel_beams';folder.mkdir(exist_ok=True)
    output=folder/'model.xml';ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert np.max(abs(data.efc_pos),initial=0)<1e-8
    for row in rows:
        for index in (0,row['specification']['segments']):
            site=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"{row['head']}__node_{index}")
            expected=np.array(row['path_source_mm'][0 if index==0 else -1])*.001
            assert np.allclose(data.site_xpos[site],expected,atol=1e-9)
    save(folder/'model_check.json',dict(parent=provenance(parent_file),assembly=provenance(assembly_file),
        model=provenance(output),generator=provenance(ROOT/'tools/build_channel_first_elbow_channel_beams.py'),
        beam_generator=provenance(ROOT/'tools/simulate_channel_first_elastic_conduits.py'),channels=rows,
        nq=model.nq,nv=model.nv,neutral_equality_max_error_m=float(np.max(abs(data.efc_pos),initial=0)),
        scope='Curved elastic-channel construction only, not driven acceptance. End clamp/translation boundary assumptions require local interface review. No manufacturing shape replacement. No guide/bone contact, self-contact, cable friction, full lattice or cross-section collapse; bone/mount geometry preserved. Other25 tendons remain endpoint surrogates.'))
    print('Elastic elbow channels:',model.nv,'velocity DOFs; neutral end constraints verified',flush=True)


if __name__=='__main__':main()
