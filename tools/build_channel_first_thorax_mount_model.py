"""Place the sternum reaction frame in the common model, without inventing supply forces."""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    previous_file=base/'elbow_supply_interface/model_check.json';previous=read(previous_file)
    old=read(ROOT/previous['assembly']['file'])
    assembly_file=base/'thorax_reaction_mount/assembly.json';assembly=read(assembly_file)
    model_file=ROOT/previous['model']['file'];assert provenance(model_file)==previous['model']
    root=ET.parse(model_file).getroot();asset=root.find('asset');thorax=root.find(".//body[@name='thorax']")
    assert thorax is not None
    owners=list(previous['geom_owners'])
    for key,source in assembly['sources'].items():
        if key in old['sources'] and source==old['sources'][key]:continue
        position=np.array(source['translation_source_mm'])*.001
        if key in old['sources']:
            row=next(r for r in owners if r['source_key']==key);assert row['owner']=='thorax'
            geom=root.find(f".//geom[@name='{row['geom']}']");geom.set('pos',vector(position))
            root.find(f"./asset/mesh[@name='{geom.get('mesh')}']").set('file',str(ROOT/source['file']))
        else:
            name='thorax_reaction_'+str(len(owners))
            ET.SubElement(asset,'mesh',name=name,file=str(ROOT/source['file']),scale='.001 .001 .001')
            ET.SubElement(thorax,'geom',name=name,type='mesh',mesh=name,pos=vector(position),
                rgba='.15 .55 .9 1',contype='0',conaffinity='0')
            owners.append(dict(source_key=key,owner='thorax',geom=name))
    for i,(key,source) in enumerate(assembly['nominal_conduit_stub_sources'].items()):
        name=f'thorax_reference_stub_{i}';ET.SubElement(asset,'mesh',name=name,file=str(ROOT/source['file']),scale='.001 .001 .001')
        ET.SubElement(thorax,'geom',name=name,mesh=name,type='mesh',mass='0',
            pos=vector(np.array(source['translation_source_mm'])*.001),rgba='.1 .8 .2 .5',contype='0',conaffinity='0')
    folder=assembly_file.parent;output=folder/'model.xml'
    ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert model.nv==previous['nv'] and len(owners)==len(assembly['sources'])==235
    save(folder/'model_check.json',dict(parent=previous['parent'],geometry_parent=provenance(previous_file),
        model=provenance(output),assembly=provenance(assembly_file),geom_owners=owners,channels=previous['channels'],
        nq=model.nq,nv=model.nv,geoms=model.ngeom,
        generator=provenance(ROOT/'tools/build_channel_first_thorax_mount_model.py'),
        scope='Fixed-thorax hardware geometry only; no full supply cable/sheath or motor. Does not add stiff or world-fixed cable via points. Thorax fixed to world remains an idealization, base load path absent. Original14 anatomical junction overlaps remain explicitly recorded; not whole chest mechanical acceptance. Local elbow beam physics inherited, not a new external supply drive result.'))
    print('Thorax port geometry model',model.nv,'DOFs',model.ngeom,'geoms;235rigidowners',flush=True)


if __name__=='__main__':main()
