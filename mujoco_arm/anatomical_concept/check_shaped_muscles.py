"""Exported geometry and nominal wire clearance; no deformability/drive pass."""
from parameterize_elbow_muscles import *
from build_shaped_elbow_muscles import capsule

def run(muscles=('FJ1486','FJ1477')):
    for fid in muscles:
        folder=DEST/fid;data=json.loads((folder/'cad_input.json').read_text(encoding='utf8'));design=json.loads((folder/'design.json').read_text(encoding='utf8'));rows=[];solids={}
        for part in data['parts']:
            name=part['name'];v,f=read_obj(str(folder/(name+'.obj')));s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
            solids[name]=s
            raw=(folder/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);stl=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)))
            volumes=[float(a.volume()) for a in s.decompose()];error=abs(stl.volume()-s.volume())/max(1,abs(s.volume()))
            rows.append({'part':name,'OBJ_status':str(s.status()),'STL_status':str(stl.status()),'decomposed_signed_volumes_mm3':volumes,'positive_material_components':sum(a>1e-6 for a in volumes),'negative_cavity_shells':sum(a< -1e-10 for a in volumes),'STL_relative_volume_error':error,'roundtrip_pass':s.status()==stl.status()==mf.Error.NoError and error<1e-4})
            if name=='TPU_muscle' or not rows[-1]['roundtrip_pass']:print(fid,name,'volumes',volumes[:10],'roundtrip',rows[-1]['roundtrip_pass'],flush=True)
        nodes=np.array(design['wire_polyline_display_mm']);wire=mf.Manifold.batch_boolean([capsule(a,b,.3) for a,b in zip(nodes,nodes[1:])],mf.OpType.Add);wire_hits=[];clearance=2.
        for name,s in solids.items():
            overlap=max(0,float((s^wire).volume()))
            if overlap>1e-6:wire_hits.append({'part':name,'volume_mm3':overlap})
            clearance=min(clearance,float(s.min_gap(wire,2)))
        result={'scope':'manufacturing geometry topology and OBJ/STL roundtrip only, NOT force, printing or actuation acceptance','parts':rows,'all_roundtrips_pass':all(r['roundtrip_pass'] for r in rows),'single_connected_material_each':all(r['positive_material_components']==1 for r in rows),'frozen':False}
        result.update(nominal_steel_diameter_mm=.6,nominal_steel_interferences=wire_hits,nominal_steel_min_clearance_mm_capped_2=clearance,nominal_steel_threading_geometry_pass=not wire_hits,neutral_part_overlaps=design['neutral_module_overlaps'])
        (folder/'geometry_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8')
        print(fid,'parts',len(rows),'steel clear',not wire_hits,'min gap mm',clearance,flush=True)

if __name__=='__main__':run()
