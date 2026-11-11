"""Explicit annular rest shells on selected seat-to-seat paths; not flexible validation."""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,export,mesh_geom
from build_channel_first_inventory import provenance
from build_channel_first_body_chain import vector


def swept_wall(points,outer=1.5,inner=.6,sides=32):
    p=np.array(points);edges=np.diff(p,axis=0);unit=edges/np.linalg.norm(edges,axis=1)[:,None]
    tangents=np.vstack([unit[0],unit[:-1]+unit[1:],unit[-1]])
    tangents/=np.linalg.norm(tangents,axis=1)[:,None]
    first=np.eye(3)[int(np.argmin(abs(tangents[0])))];previous=first-tangents[0]*(first@tangents[0]);previous/=np.linalg.norm(previous)
    bases=[]
    for t in tangents:
        radial=previous-t*(previous@t);radial/=np.linalg.norm(radial)
        bases.append((radial,np.cross(t,radial)));previous=radial
    theta=np.linspace(0.,2*np.pi,sides,endpoint=False);vertices=[]
    for radius in (outer,inner):
        for centre,(u,v) in zip(p,bases):
            vertices.extend(centre+radius*(np.cos(theta)[:,None]*u+np.sin(theta)[:,None]*v))
    faces=[];n=len(p)*sides
    for offset,reverse in [(0,False),(n,True)]:
        for i in range(len(p)-1):
            for j in range(sides):
                k=(j+1)%sides;a=offset+i*sides+j;b=offset+i*sides+k;c=offset+(i+1)*sides+k;d=offset+(i+1)*sides+j
                faces.extend([(a,c,b),(a,d,c)] if reverse else [(a,b,c),(a,c,d)])
    for j in range(sides):
        k=(j+1)%sides
        faces.extend([(j,n+j,n+k),(j,n+k,k)])
        a=(len(p)-1)*sides+j;b=(len(p)-1)*sides+k
        faces.extend([(a,b,n+b),(a,n+b,n+a)])
    vertices=np.array(vertices);faces=np.array(faces,dtype=np.uint32)
    _,boundary,nonmanifold=mesh_geom.edge_census(faces);assert boundary==nonmanifold==0
    solid=A.I._manifold().Manifold(A.I._manifold().Mesh64(vertices,faces))
    assert str(solid.status()).endswith('NoError') and solid.volume()>0
    return solid


def main():
    import mujoco
    sys.path.insert(0,str(ROOT));from project_runtime import load_mujoco_model
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    screen_file=base/'thorax_elbow_supply/route_screen.json';screen=read(screen_file)
    assembly_file=ROOT/screen['assembly']['file'];assert provenance(assembly_file)==screen['assembly']
    assembly=read(assembly_file);folder=screen_file.parent
    sources={};records=[];solids={k:load_component(s) for k,s in assembly['sources'].items()}
    bounds=np.array([s.bounding_box() for s in solids.values()]);keys=list(solids)
    for row in screen['routes']:
        selected=row['selected'];assert selected is not None
        path=np.array(selected['path_source_mm']);wall=swept_wall(path)
        actual,source=export(wall,folder/(row['id']+'_reference_sheath.obj'),path.mean(axis=0))
        collisions=[];b=np.array(actual.bounding_box())
        admitted=np.all(b[:3]<=bounds[:,3:]+1e-8,axis=1)&np.all(bounds[:,:3]<=b[3:]+1e-8,axis=1)
        for index in np.flatnonzero(admitted):
            key=keys[index];v=float(A.I.intersection_volume(actual,solids[key]))
            if v>.001:collisions.append(dict(target=key,overlap_mm3=v))
        assert not collisions,collisions
        sources['reference_sheath:'+row['id']]=source
        records.append(dict(id=row['id'],source=source,path_length_mm=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()),
            outer_radius_mm=1.5,inner_radius_mm=.6,circumferential_sides=32,centreline_rings=len(path),
            neutral_wall_interferences=collisions,
            scope='Annular reference shell with flat physical seat faces. Not a measured material, compression-resistant sheath, manufactured final muscle or validated bending result.'))
    pair_volume=float(A.I.intersection_volume(load_component(records[0]['source']),load_component(records[1]['source'])))
    assert pair_volume<=.001
    save(folder/'sheath_geometry.json',dict(route_screen=provenance(screen_file),assembly=provenance(assembly_file),
        generator=provenance(ROOT/'tools/build_channel_first_external_sheaths.py'),sources=sources,records=records,
        shell_pair_overlap_mm3=pair_volume,status='neutral_hollow_rest_geometry_not_dynamic_or_material_approval'))
    # A separate zero-DOF scene prevents these reference curves being mistaken
    # for constrained guides in the moving294DOF drive model.
    root=ET.Element('mujoco',model='neutral_supply_geometry_reference_not_dynamics')
    ET.SubElement(root,'compiler',angle='radian');ET.SubElement(root,'option',gravity='0 0 0')
    visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1100',offheight='900')
    asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody')
    ET.SubElement(world,'light',pos='-.1 -.5 1.6',dir='0 1 -1')
    all_sources={**assembly['sources'],**assembly['display_only_TPU_sources'],**sources}
    for i,(key,source) in enumerate(all_sources.items()):
        rgba='.8 .82 .76 .6'
        if key.startswith(('thorax_mount:','thorax_port:')):rgba='.1 .55 .95 1'
        elif key.startswith('tpu:'):rgba='.85 .08 .08 1'
        elif key.startswith('reference_sheath:'):rgba='.08 .7 .18 .8'
        elif key.startswith(('elbow:','supply:','stop:')):rgba='.65 .7 .75 1'
        name=f'part_{i}';ET.SubElement(asset,'mesh',name=name,file=str(ROOT/source['file']),scale='.001 .001 .001')
        ET.SubElement(world,'geom',name=name,type='mesh',mesh=name,mass='0',
            pos=vector(np.array(source.get('translation_source_mm') or [0.,0.,0.])*.001),
            rgba=rgba,contype='0',conaffinity='0')
    output=folder/'neutral_geometry.xml';ET.ElementTree(root).write(output,encoding='utf-8',xml_declaration=True)
    model=load_mujoco_model(output);assert model.nv==0
    save(folder/'neutral_scene.json',dict(model=provenance(output),sources=all_sources,
        sheath_geometry=provenance(folder/'sheath_geometry.json'),nv=model.nv,geoms=model.ngeom,
        scope='Static zero-DOF geometry scene only: actual235rigid parts, original two elbow TPU meshes, two neutral sheath shells. No fixed routing via sites in a dynamics model. Does not demonstrate deformation, full drive, shoulder muscle bodies or attached motor.'))
    print('Two reimported hollow shells clear235rigid parts; static scene',model.ngeom,'geoms,0DOF',flush=True)


if __name__=='__main__':main()
