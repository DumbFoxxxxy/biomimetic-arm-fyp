"""Direct prototype neutral all-pair solids, actual seats and internal steel."""
from build_direct_muscle_yoke import *
from build_shaped_elbow_muscles import capsule

def run():
    a=json.loads((DEST/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};hashes={str(ROOT/r['file']):hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest() for r in rows.values()};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};hard=[];soft=[]
    for (ka,sa),(kb,sb) in itertools.combinations(world.items(),2):
        v=intersection(sa,sb)
        if v>.001:(soft if 'TPU' in [rows[ka]['kind'],rows[kb]['kind']] else hard).append({'pair':[ka,kb],'volume_mm3':v})
    yoke=world['direct_muscle_yoke'];forearm=world['FJ3391'];j=json.loads((BASE/'design.json').read_text());n=np.array(j['bone_anchors'][0]['bolt_axis']);seats=[]
    for i,anchor in enumerate(j['bone_anchors']):
        q=np.array(anchor['front_spotface_mm']);patch=rod(q-n*.2,q,4.8)-rod(q-n,q+n,1.7);existing=patch^forearm;gap=float(existing.min_gap(yoke,1)) if existing.volume()>.001 else None;closing=intersection(existing,yoke.translate(-n*.01));washer=world[f'bone_{i}_front_washer'];seats.append({'anchor':i,'existing_support_patch_mm3':float(existing.volume()),'patch_to_yoke_gap_mm':gap,'closing_0_01mm_overlap_mm3':closing,'washer_to_yoke_gap_mm':float(washer.min_gap(yoke,1)),'nominal_contact_pass':existing.volume()>.1 and gap is not None and gap<.0001 and closing>.01})
    steel=[]
    for fid in ['FJ1486','FJ1477']:
        folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());off=np.array(next(p for p in a['placements'] if p['id']==fid)['translation_display_mm']);nodes=np.array(d['wire_polyline_display_mm'])+off;wire=mf.Manifold.batch_boolean([capsule(p,q,.3) for p,q in zip(nodes,nodes[1:])],mf.OpType.Add);hits=[];minimum=2
        for k,s in world.items():
            volume=intersection(wire,s)
            if volume>.000001:hits.append({'part':k,'volume_mm3':volume})
            aa=np.array(wire.bounding_box());bb=np.array(s.bounding_box())
            if np.linalg.norm(np.maximum(0,np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])))<2:minimum=min(minimum,float(wire.min_gap(s,2)))
        steel.append({'muscle':fid,'nominal_threading_path_hits':hits,'minimum_clearance_mm_capped2':minimum,'loaded_wire_or_crimp_verified':False})
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    result={'scope':'all107 manufactured neutral instances, mating seat witnesses and two internal assembly-threading paths; no loaded TPU/steel or strength proof','neutral_rigid_interferences':hard,'neutral_TPU_review_flags':soft,'bone_mount_seats':seats,'steel_neutral_threading':steel,'source_sha256':hashes,'complete_TPU_deformation_or_drive_verified':False,'frozen':False};(DEST/'interface_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('direct neutral hard',len(hard),'TPU flags',len(soft),'seat contacts',[r['nominal_contact_pass'] for r in seats],'steel hits',[len(r['nominal_threading_path_hits']) for r in steel],flush=True)

if __name__=='__main__':run()
