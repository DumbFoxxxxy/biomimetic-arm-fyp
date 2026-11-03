"""Fit existing flexible bellies around actual proximal mounting hardware.

Keeps intended pad/support contact, relieves only hardware that actually
penetrates a belly, protects the existing wire-wall material and checks the
complete readback assembly. Does not prove clip force, preload or drive.
"""
import itertools
import json
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_shoulder_origin_mounts as B
import audit_shoulder_retained_joint_fit as A
from build_shoulder_unified_muscle import local_export
from fit_shoulder_deltoid_to_collar import reconnect_outside_interface

DEST=B.DEST/'fit_revision'


def reroute_cut_segments(solid,obstacle,bone,guards,normal_at):
    """Replace a clipped local bar with a full-R1 exterior printed arch.

Both new ends join the surviving main body on opposite sides of the old
segment's principal direction. This is a recorded local lattice redesign,
not deletion of a useful fragment to manufacture a connectivity pass.
"""
    mf=A.I._manifold();changes=[]
    for old in sorted([c for c in solid.decompose() if c.volume()>1e-9],key=lambda c:-c.volume())[1:]:
        if old.volume()>5.:raise ValueError('Cut region needs full lattice regeneration')
        original_vertices=np.array(old.to_mesh64().vert_properties)[:,:3]
        centre=original_vertices.mean(0);_,_,V=np.linalg.svd(original_vertices-centre,full_matrices=False);axis=V[0]
        main=max(solid.decompose(),key=lambda c:c.volume());pv=np.array(main.to_mesh64().vert_properties)[:,:3]
        distance=np.linalg.norm(pv-centre,axis=1);projection=(pv-centre)@axis
        sides=[]
        for sign in [-1,1]:
            candidates=np.where((distance<15.)&(projection*sign>1.5))[0]
            candidates=candidates[np.argsort(distance[candidates])];safe=[]
            for index in candidates[:160]:
                p=pv[index];probe=mf.Manifold.sphere(1.,16).translate(p)
                if A.I.intersection_volume(probe,obstacle)>.001:continue
                if A.I.intersection_volume(probe,bone)>.001:continue
                if any(A.I.intersection_volume(probe,g)>.001 for g in guards):continue
                safe.append(p)
                if len(safe)>=20:break
            sides.append(safe)
        pairs=sorted((np.linalg.norm(p-q),i,j) for i,p in enumerate(sides[0]) for j,q in enumerate(sides[1]))
        chosen=None
        for _,i,j in pairs:
            p,q=sides[0][i],sides[1][j];normal=A.D.unit(normal_at((p+q)/2))
            for shift in [0.,1.5,3.,4.5,6.]:
                path=np.array([p,p+normal*shift,q+normal*shift,q]);parts=[]
                parts.extend(B.rod(a,b,1.) for a,b in zip(path,path[1:]) if np.linalg.norm(b-a)>1e-8)
                parts.extend(mf.Manifold.sphere(1.,16).translate(x) for x in path)
                rib=mf.Manifold.batch_boolean(parts,mf.OpType.Add)
                if A.I.intersection_volume(rib,obstacle)>.001:continue
                if A.I.intersection_volume(rib,bone)>.001:continue
                if any(A.I.intersection_volume(rib,g)>.001 for g in guards):continue
                end_contacts=[A.I.intersection_volume(mf.Manifold.sphere(1.,16).translate(x),main) for x in [p,q]]
                if min(end_contacts)<.02:continue
                trial=(solid-old)+rib
                if len([c for c in trial.decompose() if c.volume()>1e-9])!=1:continue
                chosen=(trial,path,rib,end_contacts);break
            if chosen is not None:break
        if chosen is None:raise ValueError('No full-radius replacement arch; retain rejected geometry')
        before=solid;solid,path,rib,contacts=chosen
        changes.append(dict(replaced_cut_segment_mm3=float(old.volume()),old_bounds_source_mm=list(old.bounding_box()),
            radius_mm=1.,path_source_mm=path.tolist(),new_material_mm3=float((solid-before).volume()),
            node_overlap_mm3=contacts,method='local full-diameter bar reroute between two surviving nodes'))
    return solid,changes


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    source=B.DEST/'assembly_check.json';parent=json.loads(source.read_text(encoding='utf8'))
    if parent['script_sha256']!=A.hashed(Path(B.__file__).resolve()):raise ValueError('Stale origin mount generator')
    solids={};bone_tri=None
    for k,v in parent['sources'].items():
        solids[k],tri,_=A.load_checked(A.ROOT/v['file'],v['translation_source_mm'],v['sha256'])
        if k=='bone:scapula':bone_tri=tri
    sources=dict(parent['sources'])
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(routes_file.read_text(encoding='utf8'))['routes']
    guards={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),1.61) for r in routes}
    lumens={r['head']:A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.86) for r in routes}
    changes={}
    for key in [k for k in solids if k.startswith('tpu:')]:
        original=solids[key]
        hit_keys=[k for k,s in solids.items() if k.startswith('origin:') and A.I.intersection_volume(original,s)>.001]
        if not hit_keys:continue
        obstacle=A.I._manifold().Manifold.batch_boolean([solids[k] for k in hit_keys],A.I._manifold().OpType.Add)
        keepout=B.clearance(obstacle)
        revised=original-keepout
        required=[]
        for r in routes:
            if 'tpu:'+r['module']==key:
                patch=(guards[r['head']]-lumens[r['head']])^original
                revised=revised+patch;required.append(patch)
        # Isolated cut wedges cannot carry load. Only remove sub-0.01mm3
        # residues outside both attachment zones and the full R1.8 channel;
        # keep the meaningful island and give it a real printed connection.
        components=sorted([c for c in revised.decompose() if c.volume()>1e-9],key=lambda c:-c.volume())
        scraps=[c for c in components[1:] if c.volume()<.01]
        removed=[]
        own_routes=[r for r in routes if 'tpu:'+r['module']==key]
        if sum(c.volume() for c in scraps)>=.01:raise ValueError('Cut residue exceeds bounded cleanup')
        for fragment in scraps:
            centre=np.mean(np.array(fragment.to_mesh64().vert_properties)[:,:3],axis=0)
            for r in own_routes:
                if A.I.intersection_volume(fragment,A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),1.8))>1e-7:
                    raise ValueError('Small fragment belongs to protected channel')
                if any(np.linalg.norm(centre-np.array(p['inner_face_source_mm']))<14 for p in r['ports'].values()):
                    raise ValueError('Small fragment belongs to attachment zone')
            revised=revised-fragment
            removed.append(dict(volume_mm3=float(fragment.volume()),bounds_source_mm=list(fragment.bounding_box()),
                reason='isolated cut wedge outside channel and both attachment zones'))
        # Allow intended contact surfaces already present; the rib must stay
        # outside actual bone/hardware and all protected wire-wall zones.
        rerouted=[]
        try:
            all_mounts=A.I._manifold().Manifold.batch_boolean([s for k,s in solids.items() if k.startswith('origin:')],A.I._manifold().OpType.Add)
            revised,ribs=reconnect_outside_interface(revised,all_mounts,solids['bone:scapula'],list(guards.values()),
                detour_direction=lambda p:A.mesh_geom.signed_distance_and_normal(p,bone_tri)[1])
        except ValueError as initial_error:
            try:
                revised,rerouted=reroute_cut_segments(revised,all_mounts,solids['bone:scapula'],list(guards.values()),
                    lambda p:A.mesh_geom.signed_distance_and_normal(p,bone_tri)[1])
                ribs=[]
            except ValueError as error:
                rejected=dict(scope='local proximal hardware relief rejected before acceptance',
                    reason=str(error),initial_reconnection_reason=str(initial_error),
                    source_sha256=A.hashed(source),script_sha256=A.hashed(Path(__file__).resolve()),
                    muscle=key,cut_obstacles=hit_keys,
                    positive_components=[dict(volume_mm3=float(c.volume()),bounds_source_mm=list(c.bounding_box()))
                                         for c in revised.decompose() if c.volume()>1e-9],stage_pass=False)
                (DEST/('rejected_'+key.split(':')[1]+'.json')).write_text(json.dumps(rejected,indent=2)+'\n',encoding='utf8',newline='\n')
                raise
        path=DEST/(key.split(':')[1]+'_hardware_relief.stl')
        origin=np.array(sources[key]['translation_source_mm'])
        actual,meta=local_export(revised,path,origin)
        cavity_repairs=[]
        if not meta['mesh_checks_pass'] and meta['positive_material_components']==1 and not meta['boundary_edges']:
            # Heal only a new microscopic edge-touching cavity. Protect every
            # cavity from the source design and every wire lumen; any minute
            # added material outside the source is explicitly recorded.
            original_voids=[]
            for old_shell in original.decompose():
                if old_shell.volume()>=0:continue
                ov,of=A.mesh_geom.manifold_mesh_to_arrays(old_shell.to_mesh64())
                original_voids.append(A.I._manifold().Manifold(A.I._manifold().Mesh64(
                    np.array(ov,dtype=np.float64,order='C',copy=True),
                    np.array(of[:,::-1],dtype=np.uint64,order='C',copy=True))))
            for shell in actual.decompose():
                volume=float(shell.volume())
                if not -1e-6<volume<0:continue
                cv,cf=A.mesh_geom.manifold_mesh_to_arrays(shell.to_mesh64())
                if np.linalg.norm(cv.max(0)-cv.min(0))>.2:continue
                filling=A.I._manifold().Manifold(A.I._manifold().Mesh64(
                    np.array(cv,dtype=np.float64,order='C',copy=True),np.array(cf[:,::-1],dtype=np.uint64,order='C',copy=True)))
                source_material=A.I.intersection_volume(filling,original)
                source_void_overlap=sum(A.I.intersection_volume(filling,v) for v in original_voids)
                if source_void_overlap>1e-9:continue
                if any(A.I.intersection_volume(filling,x)>1e-9 for x in lumens.values()):continue
                if A.I.intersection_volume(filling,all_mounts)>1e-9:continue
                actual=actual+filling
                cavity_repairs.append(dict(volume_mm3=abs(volume),bounds_source_mm=list(shell.bounding_box()),
                    source_material_mm3=source_material,original_designed_void_overlap_mm3=source_void_overlap,
                    reason='local sub-microscopic edge-cavity fill; original designed voids and lumens retained'))
            if cavity_repairs:
                archive=DEST/'iterations/submicroscopic_export_cavity';archive.mkdir(parents=True,exist_ok=True)
                if not (archive/path.name).exists():shutil.copyfile(path,archive/path.name)
                actual,meta=local_export(actual,path,origin)
        if not meta['mesh_checks_pass']:
            (DEST/('rejected_'+key.split(':')[1]+'.json')).write_text(json.dumps(meta,indent=2)+'\n',encoding='utf8',newline='\n')
            raise ValueError('Hardware relief breaks TPU mesh '+key)
        changes[key]=dict(cut_obstacles=hit_keys,removed_mm3=float((original-actual).volume()),
            added_mm3=float((actual-original).volume()),restoring_ribs=ribs,
            regenerated_lattice_segments=rerouted,removed_cut_wedges=removed,
            submicroscopic_cavity_repairs=cavity_repairs,
            additional_missing_protected_wall_mm3=sum(float((p-actual).volume()) for p in required),
            nominal_clearance_construction='six source-axis translations +/-0.1mm, with existing guard material restored',
            guaranteed_uniform_clearance=False)
        solids[key]=actual;sources[key]=meta
        print(key,'fit',changes[key],flush=True)
    contacts=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:contacts.append(dict(a=a,b=b,volume_mm3=volume,
            category='compliant contact review' if a.startswith('tpu:') or b.startswith('tpu:')
            or a.startswith('collar:liner') or b.startswith('collar:liner') else 'rigid/wire fit review'))
    wires=[]
    for row in routes:
        probe=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(probe,s)) for k,s in solids.items()]
        wires.append(dict(head=row['head'],intersections=hits,hits_over_threshold=[h for h in hits if h['volume_mm3']>.001]))
    report=dict(scope=__doc__,n_instances=len(sources),sources=sources,
        parent_source=source.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(source),
        script_sha256=A.hashed(Path(__file__).resolve()),routes_sha256=A.hashed(routes_file),
        proximal_mounts=parent['proximal_mounts'],bone_drilling=parent['bone_drilling'],TPU_hardware_reliefs=changes,
        intended_body_owners=parent['intended_body_owners'],free_terminal_constraints_required=parent['free_terminal_constraints_required'],
        material_intersections=contacts,actuation_wire_checks=wires,
        not_verified=parent['not_verified']+['existing fixture pads/jaws and bench fastening still incomplete',
            'surface gaps are nominal, not printed tolerance qualification'],stage_pass=False,actual_drive_range_deg=None)
    (DEST/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Fitted assembly',len(sources),'objects; overlaps >.05',
          [r for r in contacts if r['volume_mm3']>.05],
          'wire hits',sum(bool(w['hits_over_threshold']) for w in wires),flush=True)
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
