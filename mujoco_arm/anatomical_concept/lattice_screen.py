"""Screen BCC lattice occupancy inside ALL supplied anatomical muscle volumes.
Geometric PCA defines grid coordinates ONLY, never inferred muscle fibers.
No wire channel or endpoint is invented before their paths are labeled.
"""
from repair_inputs import *
import argparse,time

def rod(a,b,r):
    z=(b-a)/np.linalg.norm(b-a);x=np.cross(z,[1,0,0] if abs(z[0])<.8 else [0,1,0]);x/=np.linalg.norm(x);y=np.cross(z,x)
    return mf.Manifold.cylinder(float(np.linalg.norm(b-a))+.2,r,r,10).transform(np.c_[np.stack([x,y,z],1),a-z*.1])

def write_stl(path,s):
    m=s.to_mesh64();v=np.asarray(m.vert_properties)[:,:3];f=np.asarray(m.tri_verts);tri=v[f];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.maximum(np.linalg.norm(norm,axis=1),1e-24)[:,None];dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);arr=np.zeros(len(f),dt);arr['n']=norm;arr['v']=tri;path.write_bytes(b'ANATOMICAL VOLUME CLIPPED LATTICE SCREEN; NOT DRIVE READY'.ljust(80,b' ')+struct.pack('<I',len(f))+arr.tobytes())
    lines=['# BodyParts3D-derived geometry screening; CC BY-SA2.1 Japan; not print/drive acceptance.']+['v '+' '.join(format(float(c),'.17g') for c in p) for p in v]+['f '+' '.join(str(int(i)+1) for i in face) for face in f]
    path.with_suffix('.obj').write_text('\n'.join(lines)+'\n',encoding='utf8')

def run():
    parser=argparse.ArgumentParser();parser.add_argument('--only',nargs='*');parser.add_argument('--surface-cage',action='store_true');args=parser.parse_args();dest=OUT/('surface_cage_screen' if args.surface_cage else 'lattice_screen');dest.mkdir(exist_ok=True);inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));rows=[]
    config_path=Path(__file__).parent/'lattice_screen_parameters.json';config=json.loads(config_path.read_text(encoding='utf8'));spacing=np.array(config['cell_spacing_mm']);radius=config['strut_diameter_mm']/2
    existing=dest/'screen_report.json'
    if args.only and existing.exists():rows=json.loads(existing.read_text(encoding='utf8'))['muscles']
    for rec in inv['models']:
        if rec['kind']!='muscle' or (args.only and rec['id'] not in args.only):continue
        start=time.monotonic();v,f=read_obj(str(OUT/'repaired'/(rec['id']+'.obj')));v,f=weld_stl(v[f]);centre=v.mean(0);_,_,vt=np.linalg.svd(v-centre,full_matrices=False);z=vt[0];z*=1 if z[2]>=0 else -1;x=vt[1];y=np.cross(z,x);F=np.stack([x,y,z],1);local=(v-centre)@F
        body=mf.Manifold(mf.Mesh64(local.astype(np.float64),f.astype(np.uint64)))
        if body.status()!=mf.Error.NoError:raise RuntimeError('Invalid input '+rec['id'])
        lo=np.floor(local.min(0)/spacing).astype(int);hi=np.ceil(local.max(0)/spacing).astype(int);rods=[]
        for ix,iy,iz in itertools.product(*(range(a,b) for a,b in zip(lo,hi))):
            a=np.array([ix,iy,iz])*spacing;mid=a+spacing/2
            for corner in itertools.product([0,1],repeat=3):rods.append(rod(mid,a+np.array(corner)*spacing,radius))
        surface_edges=0
        if args.surface_cage:
            coarse=body.simplify(.3).to_mesh64();cv=np.asarray(coarse.vert_properties)[:,:3];cf=np.asarray(coarse.tri_verts);edges={tuple(sorted((int(a),int(b)))) for face in cf for a,b in zip(face,np.roll(face,-1))};surface_edges=len(edges)
            rods.extend(rod(cv[a],cv[b],radius) for a,b in edges if np.linalg.norm(cv[b]-cv[a])>.001)
        lattice=mf.Manifold.batch_boolean(rods,mf.OpType.Add)^body
        components=[s for s in lattice.decompose() if s.volume()>1e-6];volumes=sorted([float(s.volume()) for s in components],reverse=True);material=float(lattice.volume());source=float(body.volume());ratio=volumes[0]/material if material>0 and volumes else 0
        raw=lattice.transform(np.c_[F,centre-SHIFT]);write_stl(dest/(rec['id']+'.stl'),raw)
        # Explicit test: clipping does not create material outside anatomy.
        outside=max(0,float((lattice-body).volume()));row={'id':rec['id'],'name':rec['name'],'cell_spacing_mm':spacing.tolist(),'strut_diameter_mm':radius*2,'grid_orientation':'geometric PCA only, NOT fiber orientation','grid_origin_source_mm':centre.tolist(),'grid_basis_source':F.tolist(),'source_volume_mm3':source,'lattice_volume_mm3':material,'solid_fraction':material/source,'positive_components':len(components),'largest_component_fraction':ratio,'component_volumes_mm3':volumes,'outside_source_volume_mm3':outside,'all_components_retained':True,'wire_channel':'not designed','outer_pattern':'not designed','attachments':'not designed','print_or_drive_ready':False,'elapsed_s':time.monotonic()-start}
        row['surface_edge_struts']=surface_edges;row['outer_pattern']='triangular surface cage from .3mm simplified anatomy mesh, clipped inside exact source; exploratory, not NURBS pattern replication' if args.surface_cage else 'not designed'
        rows=[r for r in rows if r['id']!=row['id']]+[row];report={'scope':'geometry-only grid/volume screening across supplied anatomy, not full paper reproduction or contraction simulation','parameters':config,'muscles':rows};existing.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(rec['id'],'components',len(components),'largest',round(ratio,3),'fraction',round(material/source,3),flush=True)
    tree=ET.parse(OUT/'reference.xml');r=tree.getroot();r.set('model','anatomical_lattice_occupancy_screen');r.find('compiler').set('meshdir','.')
    ids={row['id'] for row in rows}
    for mesh in r.findall('./asset/mesh'):mesh.set('file',mesh.get('file') if mesh.get('name') in ids else '../assets/'+mesh.get('file'))
    for geom in r.findall('./worldbody/geom'):
        if geom.get('name') in ids:geom.set('rgba','.8 .12 .1 1')
    ET.indent(r);tree.write(dest/'reference.xml',encoding='utf-8',xml_declaration=True)

if __name__=='__main__':run()
