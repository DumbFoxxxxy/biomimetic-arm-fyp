"""Local static TPU relief for an actual shaft collar comparison.

Preserves existing cable wall guards and requires a single positive material
component after STL readback. No muscle terminal is removed or reclassified
as a tendon. The collar has not yet replaced the two old distal interfaces.
"""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_shoulder_shaft_collar as B
import audit_shoulder_retained_joint_fit as A
from build_shoulder_unified_muscle import local_export
from scipy.spatial import cKDTree


def reconnect_outside_interface(solid,keepout,bone,guards,detour_direction=None):
    """Add real R1 lattice ribs; never drop a disconnected useful component."""
    mf=A.I._manifold();repairs=[]
    components=sorted([c for c in solid.decompose() if c.volume()>1e-9],key=lambda c:-c.volume())
    if len(components)>4:raise ValueError('Relief needs a larger muscle redesign')
    primary=components[0]
    for fragment in components[1:]:
        pv=np.array(primary.to_mesh64().vert_properties)[:,:3]
        fv=np.array(fragment.to_mesh64().vert_properties)[:,:3]
        # Nearest cut faces necessarily touch the keepout; searching only
        # those endpoints misses free lattice ends outside the collar. First
        # require endpoint clearance for the actual full-radius added rib.
        def endpoint_free(p):
            probe=mf.Manifold.sphere(1.,16).translate(p)
            return (A.I.intersection_volume(probe,keepout)<=.001
                    and A.I.intersection_volume(probe,bone)<=.001
                    and all(A.I.intersection_volume(probe,g)<=.001 for g in guards))
        safe_fv=[p for p in fv if endpoint_free(p)]
        if not safe_fv:raise ValueError('No full-radius exterior endpoint on the fragment')
        fv=np.array(safe_fv)
        distances,indices=cKDTree(pv).query(fv,k=min(64,len(pv)))
        pool=sorted((float(distances[i,j]),i,int(indices[i,j]))
                    for i in range(len(fv)) for j in range(distances.shape[1]))
        q_free={};candidates=[]
        for distance,i,j in pool:
            if distance>12:break
            if j not in q_free:q_free[j]=endpoint_free(pv[j])
            if q_free[j]:candidates.append((distance,i,j))
            if len(candidates)>=80:break
        chosen=None
        centre=A.frames.joint_centres('source')['GH'];axis=A.frames.gh_probe_axes()['longitudinal']
        for length,i,j in candidates[:80]:
            if length>12:break
            p,q=fv[i],pv[j]
            midpoint=(p+q)/2;radial=midpoint-centre-axis*((midpoint-centre)@axis)
            if detour_direction is not None:radial=detour_direction(midpoint)
            radial=A.D.unit(radial)
            for shift in [0.,1.5,3.,4.5]:
                path=np.array([p,p+radial*shift,q+radial*shift,q])
                pieces=[B.rod(a,b,1.) for a,b in zip(path,path[1:]) if np.linalg.norm(b-a)>1e-8]
                pieces.extend(mf.Manifold.sphere(1.,16).translate(v) for v in path)
                rib=mf.Manifold.batch_boolean(pieces,mf.OpType.Add)
                if A.I.intersection_volume(rib,keepout)>.001:continue
                if A.I.intersection_volume(rib,bone)>.001:continue
                if any(A.I.intersection_volume(rib,g)>.001 for g in guards):continue
                trial=solid+rib
                if len([c for c in trial.decompose() if c.volume()>1e-9])>=len(components):continue
                chosen=(rib,path,shift);break
            if chosen is not None:break
        if chosen is None:raise ValueError('No local exterior rib found; preserve rejected geometry')
        rib,path,shift=chosen
        before=solid.volume();solid=solid+rib
        repairs.append(dict(fragment_volume_mm3=float(fragment.volume()),radius_mm=1.,
            path_source_mm=path.tolist(),radial_detour_mm=shift,added_mm3=float(solid.volume()-before)))
        primary=max(solid.decompose(),key=lambda c:c.volume())
    return solid,repairs


def run(folder):
    out=folder/'deltoid_relief';out.mkdir(parents=True,exist_ok=True)
    interface_file=folder/'interface_check.json'
    interface=json.loads(interface_file.read_text(encoding='utf8'))
    source=A.ROOT/interface['assembly_source']
    assembly=json.loads(source.read_text(encoding='utf8'))
    if A.hashed(source)!=interface['assembly_sha256']:raise ValueError('Stale assembly')
    parts={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0]
           for k,r in interface['sources'].items()}
    row=assembly['sources']['tpu:deltoid']
    original,_,_=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
    expanded=[]
    for body in parts.values():
        expanded.append(body)
        expanded.extend(body.translate(np.eye(3)[axis]*sign*.2)
                        for axis in range(3) for sign in [-1,1])
    keepout=A.I._manifold().Manifold.batch_boolean(expanded,A.I._manifold().OpType.Add)
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(routes_file.read_text(encoding='utf8'))['routes']
    guards=[];guard_solids=[]
    for route in routes:
        guard=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),1.61)
        volume=A.I.intersection_volume(keepout,guard)
        guards.append(dict(head=route['head'],keepout_guard_overlap_mm3=volume))
        guard_solids.append(guard)
    if any(g['keepout_guard_overlap_mm3']>.001 for g in guards):
        raise ValueError('Relief would enter protected channel-wall region')
    revised=original-keepout
    bone_row=assembly['sources']['bone:humerus']
    bone=A.load_checked(A.ROOT/bone_row['file'],bone_row['translation_source_mm'],bone_row['sha256'])[0]
    try:
        revised,ribs=reconnect_outside_interface(revised,keepout,bone,guard_solids)
    except ValueError as exc:
        rejected=dict(scope='rejected local deltoid relief/reconnection experiment',
            reason=str(exc),script_sha256=A.hashed(Path(__file__).resolve()),
            interface_sha256=A.hashed(interface_file),source_assembly_sha256=A.hashed(source),
            positive_components=[dict(volume_mm3=float(c.volume()),bounds_source_mm=list(c.bounding_box()))
                                 for c in revised.decompose() if c.volume()>1e-9],
            permitted_search='80 nearby vertex pairs, R1 struts, 0/1.5/3/4.5mm radial detours; keep collar/bone/guards clear',
            not_proven='other muscle layouts or lattice regeneration may work',stage_pass=False)
        (out/'rejected_connection.json').write_text(json.dumps(rejected,indent=2)+'\n',encoding='utf8',newline='\n')
        raise
    old_rejected=out/'rejected_export.json'
    if old_rejected.exists():
        archive=out/'iterations/unconnected_relief';archive.mkdir(parents=True,exist_ok=True)
        for name in ['rejected_export.json','deltoid_collar_relief.stl']:
            if (out/name).exists() and not (archive/name).exists():shutil.copyfile(out/name,archive/name)
    # Preserve the parent STL's local origin: changing it needlessly re-rounds
    # all remote lattice vertices in float32, including closed cavity shells.
    export_origin=np.array(row['translation_source_mm'] if row['translation_source_mm'] is not None else A.frames.SHIFT)
    actual,meta=local_export(revised,out/'deltoid_collar_relief.stl',export_origin)
    if not meta['mesh_checks_pass']:
        (out/'rejected_export.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf8',newline='\n')
        raise ValueError('Relieved muscle failed actual STL topology/connectedness/volume check')
    complete={}
    sources={}
    for key,r in assembly['sources'].items():
        complete[key]=actual if key=='tpu:deltoid' else A.load_checked(
            A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0]
        sources[key]=meta if key=='tpu:deltoid' else r
    complete.update({'collar:'+k:v for k,v in parts.items()})
    sources.update({'collar:'+k:v for k,v in interface['sources'].items()})
    intersections=[]
    for a,b in itertools.combinations(complete,2):
        volume=A.I.intersection_volume(complete[a],complete[b])
        if volume>1e-9:intersections.append(dict(a=a,b=b,volume_mm3=volume))
    wires=[]
    for route in routes:
        wire=A.D.outer_envelope(np.array(route['full_wire_path_source_mm']),.3)
        all_hits=[dict(part=k,volume_mm3=A.I.intersection_volume(wire,s)) for k,s in complete.items()]
        wires.append(dict(head=route['head'],intersections=all_hits,
                          hits_over_threshold=[h for h in all_hits if h['volume_mm3']>.001]))
    report=dict(scope=__doc__,n_instances=len(complete),sources=sources,
        assembly_source_sha256=A.hashed(source),interface_source_sha256=A.hashed(interface_file),
        routes_sha256=A.hashed(routes_file),script_sha256=A.hashed(Path(__file__).resolve()),
        relief=dict(construction='actual collar plus +/-0.2mm translations on source axes, Boolean difference',
            radial_gap_guaranteed=False,removed_mm3=float((original-actual).volume()),
            added_outside_original_mm3=float((actual-original).volume()),
            source_volume_mm3=float(original.volume()),revised_volume_mm3=float(actual.volume())),
        added_lattice_ribs=ribs,
        protected_wire_guard_radius_mm=1.61,guard_checks=guards,
        all_material_intersections=intersections,wire_checks=wires,
        not_verified=['collar not linked to shoulder muscle terminals; original distal caps retained',
            'full new muscle layout, routed short tendons, terminal fixation and new force Jacobian',
            'flexible deformation, strength, screw access and physical drive range'],
        stage_pass=False,actual_drive_range_deg=None)
    (out/'assembly_check.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Relief:',report['relief'],'topology',meta['mesh_checks_pass'],'instances',len(complete),
          'wire hits',sum(bool(w['hits_over_threshold']) for w in wires),
          'remaining overlaps >.05',[r for r in intersections if r['volume_mm3']>.05],flush=True)
    return report


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--station',type=float,default=60.)
    args=parser.parse_args();run(B.DEST/('station_'+format(args.station,'g')))
