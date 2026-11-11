"""Actual removable stop faces in the force-driven elastic-channel model."""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save,load_component
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector,uniform_geometry_inertia


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    beam_file=base/'elbow_channel_beams/model_check.json';beam=read(beam_file)
    parent_file=ROOT/beam['parent']['file'];parent=read(parent_file)
    assembly_file=base/'elbow_extension_stop/assembly.json';assembly=read(assembly_file)
    old=read(ROOT/beam['assembly']['file'])
    source_file=ROOT/beam['model']['file'];assert provenance(source_file)==beam['model']
    root=ET.parse(source_file).getroot();asset=root.find('asset')
    bodies={b.get('name'):b for b in root.findall('.//body')}
    origins={k:np.array(v) for k,v in parent['neutral_origins_m'].items()}
    owners=list(parent['geom_owners']);added=[]
    for name in ('humerus','forearm'):
        key='bone:'+name;row=next(r for r in owners if r['source_key']==key)
        geom=root.find(f".//geom[@name='{row['geom']}']");source=assembly['sources'][key]
        root.find(f"./asset/mesh[@name='{geom.get('mesh')}']").set('file',str(ROOT/source['file']))
        geom.set('pos',vector(np.array(source['translation_source_mm'])*.001-origins[name]))
        node=bodies[name].find('inertial');mass=float(node.get('mass'))
        mass+=(load_component(source).volume()-load_component(old['sources'][key]).volume())*1e-9*1200.
        com,tensor,_=uniform_geometry_inertia(source,origins[name],mass)
        pieces=[(mass,com,tensor)]
        if name=='humerus':
            for i,(k,s) in enumerate(assembly['sources'].items()):
                if not k.startswith('stop:'):continue
                mesh=f'actual_extension_stop_{i}';ET.SubElement(asset,'mesh',name=mesh,file=str(ROOT/s['file']),scale='.001 .001 .001')
                ET.SubElement(bodies[name],'geom',name=mesh,mesh=mesh,type='mesh',
                    pos=vector(np.array(s['translation_source_mm'])*.001-origins[name]),rgba='.6 .65 .7 1',contype='0',conaffinity='0')
                owners.append(dict(geom=mesh,source_key=k,owner=name))
                density=1200. if k=='stop:arm' else 7800.
                m=load_component(s).volume()*1e-9*density;c,t,_=uniform_geometry_inertia(s,origins[name],m)
                pieces.append((m,c,t));added.append(dict(key=k,mass_kg=m,density_kg_per_m3=density))
        mass=sum(p[0] for p in pieces);com=sum(m*c for m,c,_ in pieces)/mass;tensor=np.zeros((3,3))
        for m,c,t in pieces:
            d=c-com;tensor+=t+m*((d@d)*np.eye(3)-np.outer(d,d))
        node.set('mass',str(mass));node.set('pos',vector(com))
        node.set('fullinertia',vector(tensor[[0,1,2,0,0,1],[0,1,2,1,2,2]]))
    E=np.array(assembly['origin_source_mm']);F=np.array(assembly['frame_columns']);q=np.zeros(4)
    mujoco.mju_mat2Quat(q,F.flatten())
    for name,owner,z in [('fixed_extension_face','humerus',12.),('moving_extension_face','forearm',8.)]:
        ET.SubElement(bodies[owner],'geom',name=name,type='box',size='.002 .004 .002',
            pos=vector((E+F@np.array([28.,-14.,z]))*.001-origins[owner]),quat=vector(q),
            rgba='0 1 1 .25',contype='0',conaffinity='0',mass='0')
    contact=root.find('contact')
    if contact is None:contact=ET.SubElement(root,'contact')
    ET.SubElement(contact,'pair',name='actual_extension_stop',geom1='fixed_extension_face',geom2='moving_extension_face',
        condim='1',margin='0',gap='0',solref='-100000000 -20000',solimp='.9999 .9999 .0001')
    joint=root.find(".//joint[@name='elbow']");old_limit=dict(joint.attrib);joint.set('limited','false')
    root.find('option').set('iterations','100');root.find('option').set('tolerance','1e-10')
    folder=base/'elbow_extension_stop';model_file=folder/'model.xml'
    ET.ElementTree(root).write(model_file,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(model_file);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert model.nv==beam['nv']
    save(folder/'model_check.json',dict(parent=provenance(parent_file),beam_parent=provenance(beam_file),
        assembly=provenance(assembly_file),model=provenance(model_file),channels=beam['channels'],geom_owners=owners,
        original_elbow_joint=old_limit,elbow_software_limit_enabled=False,added_mass_assumptions=added,
        contact_law=dict(solref=[-1e8,-20000],solimp=[.9999,.9999,.0001],frictionless=True),
        nq=model.nq,nv=model.nv,generator=provenance(ROOT/'tools/build_channel_first_elbow_stop_model.py'),
        scope='Only two actual stop blocks contact-solved; all other rigid collision response remains absent and must be independently screened. No angle commands. Provisional densities/inertia, stiff soft-contact law uncalibrated. Two elbow beams; other25 muscles endpoint proxies. Not full motion or hardware acceptance.'))
    print('Actual stop contact model',model.nv,'DOFs, software elbow limit disabled',flush=True)


if __name__=='__main__':main()
