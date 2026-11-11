"""Mount real supply-end hardware; external chest/sheath forces remain absent."""
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
    check_file=base/'elbow_extension_stop/model_check.json';check=read(check_file)
    source_file=ROOT/check['model']['file'];assert provenance(source_file)==check['model']
    parent=read(ROOT/check['parent']['file'])
    assembly_file=base/'elbow_supply_interface/assembly.json';assembly=read(assembly_file)
    old=read(ROOT/check['assembly']['file']);root=ET.parse(source_file).getroot()
    bodies={b.get('name'):b for b in root.findall('.//body')};asset=root.find('asset')
    origin=np.array(parent['neutral_origins_m']['humerus']);body=bodies['humerus'];node=body.find('inertial')
    owners=list(check['geom_owners'])
    mass=float(node.get('mass'));centre=np.fromstring(node.get('pos'),sep=' ')
    v=np.fromstring(node.get('fullinertia'),sep=' ');tensor=np.array([[v[0],v[3],v[4]],[v[3],v[1],v[5]],[v[4],v[5],v[2]]])
    shift=lambda c:(c@c)*np.eye(3)-np.outer(c,c)
    first=mass*centre;at_origin=tensor+mass*shift(centre);increments=[]
    def contribution(key,source,sign):
        nonlocal mass,first,at_origin
        density=1200. if key.endswith(('rigid_cap','retainer')) else 7800.
        m=load_component(source).volume()*1e-9*density
        c,t,_=uniform_geometry_inertia(source,origin,m)
        mass+=sign*m;first+=sign*m*c;at_origin+=sign*(t+m*shift(c))
        increments.append(dict(key=key,sign=sign,mass_kg=m,assumed_density_kg_per_m3=density))
    for key,source in assembly['sources'].items():
        if key in old['sources'] and source==old['sources'][key]:continue
        if key in old['sources']:
            row=next(r for r in owners if r['source_key']==key);assert row['owner']=='humerus'
            geom=root.find(f".//geom[@name='{row['geom']}']")
            root.find(f"./asset/mesh[@name='{geom.get('mesh')}']").set('file',str(ROOT/source['file']))
            geom.set('pos',vector(np.array(source['translation_source_mm'])*.001-origin))
            contribution(key,old['sources'][key],-1.)
        else:
            name='supply_hardware_'+str(len(owners));ET.SubElement(asset,'mesh',name=name,file=str(ROOT/source['file']),scale='.001 .001 .001')
            ET.SubElement(body,'geom',name=name,mesh=name,type='mesh',
                pos=vector(np.array(source['translation_source_mm'])*.001-origin),rgba='.95 .55 .08 1',contype='0',conaffinity='0')
            owners.append(dict(source_key=key,owner='humerus',geom=name))
        contribution(key,source,1.)
    centre=first/mass;tensor=at_origin-mass*shift(centre);assert np.linalg.eigvalsh(tensor).min()>0
    node.set('mass',str(mass));node.set('pos',vector(centre));node.set('fullinertia',vector(tensor[[0,1,2,0,0,1],[0,1,2,1,2,2]]))
    # Reference stubs are display-only; they are not rigid route constraints.
    for i,(key,source) in enumerate(assembly['conduit_stub_sources'].items()):
        name='reference_sheath_stub_'+str(i)
        ET.SubElement(asset,'mesh',name=name,file=str(ROOT/source['file']),scale='.001 .001 .001')
        ET.SubElement(body,'geom',name=name,mesh=name,type='mesh',mass='0',
            pos=vector(np.array(source['translation_source_mm'])*.001-origin),rgba='.1 .8 .2 .5',contype='0',conaffinity='0')
    folder=base/'elbow_supply_interface';output=folder/'model.xml'
    ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert model.nv==check['nv'] and len(owners)==208
    save(folder/'model_check.json',dict(parent=check['parent'],mechanical_parent=provenance(check_file),
        model=provenance(output),assembly=provenance(assembly_file),geom_owners=owners,channels=check['channels'],
        nq=model.nq,nv=model.nv,geoms=model.ngeom,humerus_mass_kg=mass,mass_increments=increments,
        display_only_stubs=assembly['conduit_stub_sources'],
        generator=provenance(ROOT/'tools/build_channel_first_elbow_supply_model.py'),
        scope='Geometry/mass update of local elbow beam drive. New hardware follows humerus; two sheath stubs are neutral reference displays, not force-transmitting or fixed via-point constraints. Full chest supply, sheath compression/bending and ferrule retention absent. No rigid contact response except stop faces. Base inertia remains provisional; cap signed material increments use assumed density.'))
    print('Supply hardware model',model.nv,'DOFs',model.ngeom,'geoms;208rigid owners; estimated humerus mass',mass,flush=True)


if __name__=='__main__':main()
