"""Full manufactured lattice under a prescribed geometric mapping hypothesis.

Uses original150k-vertex bodies, not decorative muscle envelopes. Cross-section
area compensation is a geometric volume-budget choice, NOT calibrated TPU
elasticity. Inverted mappings are rejected before interpreting intersections.
"""
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation,Slerp
from pathlib import Path

def run():
    assert EXTENSOR_REANCHOR
    parent=DEST/'paired_swivel_support';candidate=parent/'captive_keepers';source=candidate/'wire_pose_study/wire_guides_quick.json';study=json.loads(source.read_text());a=json.loads((candidate/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};motion=json.loads((parent/'conditional_motion.json').read_text())['poses'];dest=candidate/'full_mesh_mapping_probe';dest.mkdir(exist_ok=True);world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if r['kind']!='TPU'};inputs=[source,candidate/'assembly.json',parent/'conditional_motion.json']+[ROOT/r['file'] for r in rows.values()];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};results=[]
    angles=[0,90,150]
    for fid in ['FJ1486','FJ1477']:
        folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());hashes[str(folder/'design.json')]=hashlib.sha256((folder/'design.json').read_bytes()).hexdigest();r=rows[fid+'_TPU_muscle'];solid,V,faces=load(ROOT/r['file']);V=V+np.array(r['translation_display_mm']);offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));A,B=np.array(d['pad_centres_display_mm'])+offset;z=np.array(d['basis_columns_display'])[:,2];centres=np.array([g['centre_display_mm'] for g in d['actual_short_tube_guides']])+offset;old=np.vstack([A+z*2,centres,B-z*2]);heights=old@z;assert np.all(np.diff(heights)>0);edges=np.unique(np.sort(np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);rest_edge=np.linalg.norm(V[edges[:,1]]-V[edges[:,0]],axis=1);long=rest_edge>.05;rest_volume=float(solid.volume())
        for angle in angles:
            wire=next(p for p in study['poses'] if p['muscle']==fid and p['elbow_deg']==angle);pose=next(p for p in motion if p['elbow_deg']==angle);Tc=np.array(pose['cap_poses'][fid]);Tf=np.array(pose['forearm_pose']);new=np.vstack([A+z*2,np.array(wire['centreline_display_mm'])[2:-2],(B-z*2)@Tc[:,:3].T+Tc[:,3]]);R=np.array([np.eye(3),*[np.array(p['pose_affine_mm'])[:,:3] for p in wire['guide_checks']],Tc[:,:3]]);slerp=Slerp(heights,Rotation.from_matrix(R));stretch=np.linalg.norm(np.diff(new,axis=0),axis=1)/np.linalg.norm(np.diff(old,axis=0),axis=1);alpha=np.r_[1,1/np.sqrt(np.maximum((stretch[:-1]+stretch[1:])/2,.01)),1]
            def warp(points):
                h=points@z;t=np.clip(h,heights[0],heights[-1]);C0=np.column_stack([np.interp(t,heights,old[:,i]) for i in range(3)]);C1=np.column_stack([np.interp(t,heights,new[:,i]) for i in range(3)]);scale=np.interp(t,heights,alpha);mapped=C1+slerp(t).apply((points-C0)*scale[:,None]);mapped[h<=heights[0]]=points[h<=heights[0]];end=h>=heights[-1];mapped[end]=points[end]@Tc[:,:3].T+Tc[:,3];return mapped
            moved=warp(V);jac=np.empty((len(V),3,3));epsilon=.01
            for k in range(3):delta=np.eye(3)[k]*epsilon;jac[:,:,k]=(warp(V+delta)-warp(V-delta))/(2*epsilon)
            determinants=np.linalg.det(jac);length=np.linalg.norm(moved[edges[:,1]]-moved[edges[:,0]],axis=1);ratios=length[long]/rest_edge[long];inverted=bool(np.any(determinants<=0));record={'muscle':fid,'elbow_deg':angle,'source_vertex_count':len(V),'source_triangle_count':len(faces),'mapping_scope':'prescribed centreline/frame skinning, radius compensation1/sqrt(segmentstretch), dense2mm endpads moved rigidly with actual cap frames; NO calibrated material, force or contact solution','minimum_sampled_surface_vertex_Jacobian':float(determinants.min()),'maximum_sampled_surface_vertex_Jacobian':float(determinants.max()),'nonpositive_Jacobian_vertex_count':int(np.sum(determinants<=0)),'surface_edge_stretch_min_max':[float(ratios.min()),float(ratios.max())],'surface_edge_stretch_p01_p99':np.percentile(ratios,[1,99]).tolist(),'input_wire_stage_has_interferences':bool(wire['steel_rigid_intersections'] or any(g['rigid_intersections'] for g in wire['guide_checks'])),'TPU_contact_classification':'soft contacts require material/contact solution; neither automatic hard failure nor verified pass','actual_TPU_deformation_verified':False,'physical_drive_verified':False}
            if not inverted:
                body=mf.Manifold(mf.Mesh64(moved.astype(np.float64),faces.astype(np.uint64)));record['mapped_mesh_status']=str(body.status());record['material_volume_change_fraction']=(float(body.volume())-rest_volume)/rest_volume;hits=[]
                if body.status()==mf.Error.NoError:
                    for key,other in world.items():
                        owner=rows[key]['owner'];T=np.array(pose['cap_poses'][owner[:6]]) if owner.endswith('_distal_pivot') else Tf if owner=='elbow' else None;target=other.transform(T) if T is not None else other;volume=intersection(body,target)
                        if volume>.001:hits.append({'part':key,'volume_mm3':volume,'classification':'prescribed_TPU_contact_review'})
                    record['steel_intersection_with_mapped_TPU_mm3']=intersection(body,capsule_chain_for_probe(np.array(wire['centreline_display_mm'])));record['rigid_contact_reviews']=hits
            else:record['contact_check_skipped_reason']='mapping has local inversion; cannot interpret the warped mesh as a physical body'
            results.append(record);print(fid,'q',angle,'minJ',record['minimum_sampled_surface_vertex_Jacobian'],'invertedverts',record['nonpositive_Jacobian_vertex_count'],'contacts',len(record.get('rigid_contact_reviews',[])),flush=True);(dest/'partial_checkpoint.json').write_text(json.dumps({'scope':'prescribed fullmesh diagnostic only','poses':results},indent=2),encoding='utf8')
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());(dest/'mapping_review.json').write_text(json.dumps({'scope':'diagnostic of a prescribed whole-manufacturing-mesh mapping; passing samples would still not prove TPU equilibrium/self-contact/incompressibility/friction or actual driverange; negative results reject this mapping, not the muscle concept','source_sha256':hashes,'poses':results,'frozen':False},indent=2),encoding='utf8')

def capsule_chain_for_probe(points):
    return mf.Manifold.batch_boolean([rod(p,q,.3,24) for p,q in zip(points,points[1:])]+[mf.Manifold.sphere(.3,16).translate(p) for p in points[1:-1]],mf.OpType.Add)

if __name__=='__main__':run()
