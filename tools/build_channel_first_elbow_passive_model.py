"""Restore real cap hinges and bore load sites; no prescribed swivel trajectories."""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector,uniform_geometry_inertia


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    assembly_file=base/'elbow_mount_candidate/assembly.json';assembly=read(assembly_file)
    pad_file=base/'unified_SC_AC/scapular_pad_dynamic_screen.json';pad=read(pad_file)
    model_file=ROOT/pad['model']['file'];assert provenance(model_file)==pad['model']
    endpoint=read(ROOT/pad['endpoint_report']['file']);chain=read(ROOT/endpoint['chain_check']['file'])
    root=ET.parse(model_file).getroot();asset=root.find('asset')
    bodies={b.get('name'):b for b in root.findall('.//body')}
    origins={k:np.array(v) for k,v in chain['neutral_body_origins_m'].items()};origins['thorax']=np.zeros(3)
    pivot_rows=[]
    legacy=ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair'
    for head,folder in [('FJ1486','cap_swivel_flexor'),('FJ1477','cap_swivel_study')]:
        screen_file=legacy/folder/'screen.json';screen=read(screen_file)
        centre=(np.array(screen['pivot_zero_display_mm'])+[-140.,-70.,1100.])*.001
        name=head+'_distal_pivot';origins[name]=centre
        body=ET.SubElement(bodies['forearm'],'body',name=name,pos=vector(centre-origins['forearm']))
        ET.SubElement(body,'joint',name=head+'_swivel',type='hinge',axis=vector(screen['trunnion_axis_zero_display']),
                      stiffness='.02',damping='.001',limited='false')
        bodies[name]=body
        mass=0.;weighted=np.zeros(3);pieces=[]
        for key,ownership in assembly['legacy_owners'].items():
            if ownership['legacy_owner']!=name:continue
            source=assembly['sources'][key]
            from build_s1_teres_origin_mount import load_component
            density=1200. if key.endswith('_cap_rotor') else 7800.
            volume=load_component(source).volume()*1e-9;m=volume*density
            com,tensor,_=uniform_geometry_inertia(source,centre,m)
            mass+=m;weighted+=m*com;pieces.append((m,com,tensor))
        com=weighted/mass;tensor=np.zeros((3,3))
        for m,c,t in pieces:
            d=c-com;tensor+=t+m*((d@d)*np.eye(3)-np.outer(d,d))
        ET.SubElement(body,'inertial',mass=str(mass),pos=vector(com),
            fullinertia=vector(tensor[[0,1,2,0,0,1],[0,1,2,1,2,2]]))
        pivot_rows.append(dict(head=head,body=name,source=provenance(screen_file),
            centre_source_mm=(centre*1000).tolist(),axis_source=screen['trunnion_axis_zero_display'],
            assumed_mass_kg=mass,stiffness_Nm_per_rad=.02,damping_Nms_per_rad=.001))
    # Replace the actual mounted humerus in the existing geometry slot.
    hum_source=assembly['sources']['bone:humerus']
    hum_geom=next(r['geom'] for r in chain['geom_owners'] if r['source_key']=='bone:humerus')
    node=root.find(f".//geom[@name='{hum_geom}']")
    root.find(f"./asset/mesh[@name='{node.get('mesh')}']").set('file',str(ROOT/hum_source['file']))
    node.set('pos',vector(np.array(hum_source['translation_source_mm'])*.001-origins['humerus']))
    com,tensor,_=uniform_geometry_inertia(hum_source,origins['humerus'],.6)
    inertial=bodies['humerus'].find('inertial');inertial.set('pos',vector(com))
    inertial.set('fullinertia',vector(tensor[[0,1,2,0,0,1],[0,1,2,1,2,2]]))
    geom_owners=list(chain['geom_owners'])
    for i,(key,ownership) in enumerate(assembly['legacy_owners'].items()):
        owner={'GH':'humerus','elbow':'forearm'}.get(ownership['legacy_owner'],ownership['legacy_owner'])
        assert owner in bodies,key
        source=assembly['sources'][key];path=ROOT/source['file'];assert provenance(path)['sha256']==source['sha256']
        name=f'elbow_interface_{i}'
        ET.SubElement(asset,'mesh',name=name,file=str(path),scale='.001 .001 .001')
        position=np.array(source.get('translation_source_mm') or [0.,0.,0.])*.001-origins[owner]
        ET.SubElement(bodies[owner],'geom',name=name,type='mesh',mesh=name,pos=vector(position),
            rgba='.65 .7 .75 1',contype='0',conaffinity='0')
        geom_owners.append(dict(geom=name,source_key=key,owner=owner))
    ports=[]
    for label,muscle in assembly['muscles'].items():
        head=muscle['anatomical_id'];wire=muscle['wire_path_source_mm']
        for role,point,owner in [('origin',wire[1],'humerus'),('distal',wire[-2],head+'_distal_pivot')]:
            site=root.find(f".//site[@name='{head}_{role}']")
            for body in bodies.values():
                if site in list(body):body.remove(site);break
            bodies[owner].append(site);site.set('pos',vector(np.array(point)*.001-origins[owner]))
            ports.append(dict(head=head,role=role,owner=owner,source_position_mm=point))
    dest=base/'elbow_passive_model';dest.mkdir(exist_ok=True)
    output=dest/'model.xml';ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    assert model.nv==12
    for port in ports:
        i=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,port['head']+'_'+port['role'])
        assert np.allclose(data.site_xpos[i]*1000,port['source_position_mm'],atol=1e-7)
    save(dest/'model_check.json',dict(assembly=provenance(assembly_file),parent_pads=provenance(pad_file),
        endpoint_report=pad['endpoint_report'],model=provenance(output),pivots=pivot_rows,ports=ports,
        geom_owners=geom_owners,neutral_origins_m={k:v.tolist() for k,v in origins.items()},
        pads=pad['pads'],nq=model.nq,nv=model.nv,geoms=model.ngeom,
        generator=provenance(ROOT/'tools/build_channel_first_elbow_passive_model.py'),
        scope='Two passive cap hinges restored. Bore-to-bore endpoint force surrogate, not the complete bent cable channel. Cap inertia uses assumed densities1200/7800kg/m3; base segment masses remain provisional inclusive estimates. No prescribed swivel angle, mechanical stop, rigid collision response, friction or manufacturing approval.'))
    print('Passive model',model.nv,'DOFs',model.ngeom,'geoms; actual bore sites verified',flush=True)


if __name__=='__main__':main()
