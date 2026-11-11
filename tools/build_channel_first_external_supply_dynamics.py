"""Free external sheath beams between real housings, steel continuous into original muscle."""
import argparse
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector
from simulate_channel_first_elastic_conduits import sample_path


def beam_tree(head,points,axial_modulus,bending_modulus,spacing):
    """Clamp at actual outlet; all intermediate beam nodes are free."""
    sampled=sample_path(np.array(points),spacing)*.001;origin=sampled[0].copy();p=sampled-origin
    lengths=np.linalg.norm(np.diff(p,axis=0),axis=1);area=np.pi*(.0015**2-.0006**2)
    moment=np.pi/4*(.0015**4-.0006**4);EA=axial_modulus*1e6*area;EI=bending_modulus*1e6*moment
    prefix=head+'_supply';root=ET.Element('body',name=prefix+'__clamp',pos=vector(origin))
    ET.SubElement(root,'site',name=prefix+'__node_0',size='.00035',rgba='0 .8 .1 1')
    parent=root;axial=[];bending=[];total_mass=0.
    for i,delta in enumerate(np.diff(p,axis=0)):
        length=lengths[i];node=ET.SubElement(parent,'body',name=f'{prefix}__end_{i}',pos=vector(delta))
        ET.SubElement(node,'joint',name=f'{prefix}__axial_{i}',type='slide',axis=vector(delta/length),
            stiffness=str(EA/length),damping=str(.5*.005/length),limited='true',range=vector([-.5*length,.5*length]))
        axial.append(EA/length)
        if i<len(lengths)-1:
            lever=.5*(length+lengths[i+1]);stiffness=EI/lever;bending.append(stiffness)
            ET.SubElement(node,'joint',name=f'{prefix}__bend_{i+1}',type='ball',stiffness=str(stiffness),damping=str(.00003*.005/lever))
        lump=.5*length+(.5*lengths[i+1] if i<len(lengths)-1 else 0.)
        mass=1200.*area*lump;total_mass+=mass;inertia=mass*(lump*lump/12.+(.0015**2+.0006**2)/4.)
        ET.SubElement(node,'inertial',pos='0 0 0',mass=str(mass),diaginertia=vector([inertia]*3))
        ET.SubElement(node,'site',name=f'{prefix}__node_{i+1}',size='.00035',rgba='0 .8 .1 1')
        parent=node
    return root,parent,dict(segments=len(lengths),source_origin_m=origin.tolist(),reference_nodes_source_mm=(sampled*1000).tolist(),
        segment_rest_lengths_m=lengths.tolist(),axial_stiffness_N_per_m=axial,bending_stiffness_Nm_per_rad=bending,
        outer_radius_mm=1.5,inner_radius_mm=.6,axial_modulus_MPa=axial_modulus,bending_modulus_MPa=bending_modulus,
        EA_N=EA,EI_Nm2=EI,assumed_density_kg_m3=1200.,moving_mass_kg=total_mass,
        discretization_length_mm=float(lengths.sum()*1000),rotational_model='isotropic ball spring EI/L, including uncalibrated twist penalty',
        inertia_model='adjacent half segment masses, isotropic nodal inertia; fixed root half mass omitted',
        damping_model='unmeasured axial0.5Ns/m and rotational3e-5Nms/rad at5mm, inversely scaled by segment length',
        numerical_axial_guard_fraction=.5)


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',default='axial15_bend15');parser.add_argument('--axial-modulus-mpa',type=float,default=15.)
    parser.add_argument('--bending-modulus-mpa',type=float,default=15.);parser.add_argument('--spacing-mm',type=float,default=5.)
    parser.add_argument('--timestep',type=float,default=.00002)
    args=parser.parse_args();assert all(c.isalnum() or c in '_-' for c in args.case)
    assert min(args.axial_modulus_mpa,args.bending_modulus_mpa,args.spacing_mm,args.timestep)>0
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    previous_file=base/'thorax_reaction_mount/model_check.json';previous=read(previous_file)
    routes_file=base/'thorax_elbow_supply/route_screen.json';routes=read(routes_file)
    ports_file=base/'elbow_supply_interface/route_endpoints.json';ports=read(ports_file)
    source_file=ROOT/previous['model']['file'];assert provenance(source_file)==previous['model']
    root=ET.parse(source_file).getroot();bodies={b.get('name'):b for b in root.findall('.//body')}
    root.find('option').set('timestep',str(args.timestep))
    parent=read(ROOT/previous['parent']['file']);origins=parent['neutral_origins_m']
    equality=root.find('equality');wire_paths=[];extra=[];clamps=[]
    # Remove rigid-looking old stub displays. They are not beam supports.
    for body in bodies.values():
        for geom in list(body.findall('geom')):
            if geom.get('name','').startswith(('reference_sheath_stub_','thorax_reference_stub_')):body.remove(geom)
    for row in routes['routes']:
        head=row['id'];selected=row['selected'];assert selected and not selected['interior_fixed_guides']
        # Flexible interval starts/ends at the actual outlets. The4mm seated
        #ends remain housing-defined straight passages, not compliant spans.
        path=np.array(selected['path_source_mm'])[1:-1]
        subtree,last,spec=beam_tree(head,path,args.axial_modulus_mpa,args.bending_modulus_mpa,args.spacing_mm)
        bodies['thorax'].append(subtree);prefix=head+'_supply';count=spec['segments']
        port=next(p for p in ports['routes'] if p['id']==head)
        end=np.array(port['physical_ferrule_outlet_source_mm'])*.001
        seat=np.array(port['physical_sheath_reaction_seat_source_mm'])*.001
        hum_origin=np.array(origins['humerus'])
        for name,point in [('outlet',end),('seat',seat)]:
            ET.SubElement(bodies['humerus'],'site',name=f'{prefix}__housing_{name}',pos=vector(point-hum_origin),size='.00035')
        ET.SubElement(last,'site',name=prefix+'__engagement_inner',pos=vector(seat-end),size='.00035')
        specifications=[(prefix+'_outlet_connect',f'{prefix}__node_{count}',prefix+'__housing_outlet'),
                        (prefix+'_tangent_connect',prefix+'__engagement_inner',prefix+'__housing_seat')]
        for name,s1,s2 in specifications:
            ET.SubElement(equality,'connect',name=name,site1=s1,site2=s2,solref='.0002 1',solimp='.9999 .99999 .00001')
        chest_seat=np.array(selected['end_supports_only'][0]['seat_source_mm'])*.001
        ET.SubElement(bodies['thorax'],'site',name=prefix+'__input_seat',pos=vector(chest_seat),size='.00035')
        internal=next(r for r in previous['channels'] if r['head']==head)
        names=[prefix+'__input_seat']+[f'{prefix}__node_{i}' for i in range(count+1)]
        names += [prefix+'__housing_seat']+[f'{head}__node_{i}' for i in range(internal['specification']['segments']+1)]
        tendon=root.find(f"./tendon/spatial[@name='{head}']")
        for item in list(tendon):tendon.remove(item)
        for name in names:ET.SubElement(tendon,'site',site=name)
        wire_paths.append(dict(id=head,site_names=names,
            input='ideal cable force boundary at actual thorax ferrule seat; no motor/spool geometry or inboard tail yet',
            path='thorax4mm passage -> free external beam -> humerus4mm passage -> rigid cap/TPU end transition proxy -> original internal beam -> passive distal cap'))
        extra.append(dict(head=prefix,specification=spec,path_source_mm=path.tolist(),region='external_supply',
            proximal_support='thorax actual outlet, tangent from real4mm engagement; ideal clamp approximation',
            distal_support='two point constraints on same moving end body at actual humerus outlet and4mm-inboard seat; axial twist free',
            intermediate_support='all nodes free; no world mocap/via sites',
            modulus_note='Independent axial/bending equivalent moduli, not a selected or tested material.'))
        clamps.append(dict(id=head,equalities=[r[0] for r in specifications],
            sheath_sites=[f'{prefix}__node_{count}',prefix+'__engagement_inner'],
            housing_sites=[prefix+'__housing_outlet',prefix+'__housing_seat']))
    folder=base/'thorax_elbow_supply_dynamic'/args.case;folder.mkdir(parents=True,exist_ok=True)
    output=folder/'model.xml';ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    eq_error=float(np.max(abs(data.efc_pos[data.efc_type==mujoco.mjtConstraint.mjCNSTR_EQUALITY]),initial=0));assert eq_error<1e-8
    for row in extra:
        spec=row['specification'];sites=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,f"{row['head']}__node_{i}") for i in range(spec['segments']+1)]
        assert np.allclose(data.site_xpos[sites]*1000,spec['reference_nodes_source_mm'],atol=1e-8)
    save(folder/'model_check.json',dict(parent=previous['parent'],geometry_parent=provenance(previous_file),
        assembly=previous['assembly'],model=provenance(output),routes=provenance(routes_file),arm_ports=provenance(ports_file),
        generator=provenance(ROOT/'tools/build_channel_first_external_supply_dynamics.py'),geom_owners=previous['geom_owners'],
        channels=[*previous['channels'],*extra],wire_paths=wire_paths,external_clamps=clamps,
        nq=model.nq,nv=model.nv,neutral_equality_error_m=eq_error,parameters=vars(args),
        scope='Two external sheath beams and continuous steel routes into original two internal beams. End boundary laws are idealized4mm housing support with retention unresolved. No intermediate fixed guides. No external wall/bone/self contact, friction, lumen collapse, actual motor, distal terminal review, full TPU lattice or strength. Nominal initial tension still inherited from old endpoint allocation, not recalculated equilibrium. Independent parameters are sensitivity cases, not proof of any material.'))
    print(args.case,model.nv,'DOFs; neutral constraints',eq_error,'continuous tendons verified',flush=True)


if __name__=='__main__':main()
