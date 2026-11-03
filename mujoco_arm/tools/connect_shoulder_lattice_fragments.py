"""Add explicit TPU transition ribs to disconnected shoulder lattice pieces.

Keeps real material pieces and checks every added rib against actual bones,
caps and neighbouring muscles. This is a static connection candidate; added
ribs alter compliance and do not validate drive, strength or installation.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from scipy.spatial import cKDTree
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import check_shoulder_assembly as C
import build_shoulder_tendon_segments as T
import interference as I
import mesh_geom
from build_shoulder_unified_muscle import local_export,cap_keepout


def vertices(solid):
    return mesh_geom.manifold_mesh_to_arrays(solid.to_mesh64())[0]


def candidate_links(a,b):
    yield np.array([a,b])
    # Local arch candidates permit a rib to pass around a cap, rather than
    # crossing it or joining through its wire lumen. Geometry stays explicit.
    directions=[D.unit(v) for v in itertools.product([-1,0,1],repeat=3) if any(v)]
    for height in [2.,4.,6.,9.,12.]:
        for direction in directions:
            path=np.array([a,a+direction*height,b+direction*height,b])
            if np.linalg.norm(np.diff(path,axis=0),axis=1).sum()<=30:yield path


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',required=True,
        choices=['infraspinatus','subscapularis','teres_major'])
    ap.add_argument('--parent-report',type=Path)
    args=ap.parse_args()
    folder=D.DEST/('unified_'+args.module+'_fit_gap')
    parent_file=args.parent_report or folder/'check.json';meta=json.loads(parent_file.read_text(encoding='utf8'))
    own=meta['parts']['TPU'];origin=np.array(own['translation_source_mm'])
    def load_row(row):
        path=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Stale part '+str(path))
        return C.load_stl(path).translate(row['translation_source_mm'])
    solid=load_row(own)
    obstacles={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()}
    for name in ['deltoid','infraspinatus','subscapularis','teres_major']:
        suffix='_outlet_socket' if name=='deltoid' else '_fit_gap'
        if name!='deltoid' and (D.DEST/('unified_'+name+'_fitted/check.json')).exists():suffix='_fitted'
        report=json.loads((D.DEST/('unified_'+name+suffix)/'check.json').read_text(encoding='utf8'))
        for key,row in report['parts'].items():
            if name==args.module and key=='TPU':continue
            obstacles[name+':'+key]=load_row(row)
    source=D.frames.PROJECT_ROOT/meta['master_routes']
    if hashlib.sha256(source.read_bytes()).hexdigest()!=meta['master_sha256']:raise ValueError('Stale routes')
    routes=[r for r in json.loads(source.read_text(encoding='utf8'))['routes'] if r['module']==args.module]
    every_route={r['head']:r for r in json.loads(source.read_text(encoding='utf8'))['routes']}
    masks=[obstacles[k] for k in D.SP.BONES]
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    for seat in seats:
        port=every_route[seat['head']]['ports'][seat['end']]
        inner=np.array(port['inner_face_source_mm'])-port['offset_from_seat_mm']
        masks.append(cap_keepout(inner,seat['normal_source'],seat['cap_thickness_mm'],seat['cap_radius_mm']))
    for route in every_route.values():
        if outlet:=route.get('rigid_origin_outlet'):
            masks.append(D.outer_envelope(np.array(outlet['path_source_mm']),outlet['outer_radius_mm']+.10))
    def bbox_overlap(a,b):
        aa=np.array(a.bounding_box());bb=np.array(b.bounding_box())
        return bool(np.all(aa[:3]<=bb[3:]) and np.all(bb[:3]<=aa[3:]))
    lumens=I._manifold().Manifold.batch_boolean([
        D.outer_envelope(np.array(r['full_wire_path_source_mm']),.86) for r in routes],I._manifold().OpType.Add)
    links=[];attempts=0;unresolved=[]
    components=sorted([p for p in solid.decompose() if p.volume()>1e-9],key=lambda p:p.volume(),reverse=True)
    original_components=[dict(volume_mm3=p.volume(),bounds_source_mm=list(p.bounding_box())) for p in components]
    # Target each real initial fragment, accepting only actual unions that
    # reduce the number of positive material components.
    for fragment in components[1:]:
        now=sorted([p for p in solid.decompose() if p.volume()>1e-9],key=lambda p:p.volume(),reverse=True)
        if len(now)==1:break
        main=now[0];target=vertices(fragment);base=vertices(main)
        main_material=solid^main
        if I.intersection_volume(main_material,fragment)>.5*fragment.volume():continue
        distance,index=cKDTree(base).query(target)
        order=np.argsort(distance)
        pairs=[]
        for i in order:
            a=target[i];b=base[index[i]]
            if all(np.linalg.norm(a-x)>1. for x,_ in pairs):pairs.append((a,b))
            if len(pairs)>=12:break
        found=None
        candidates=[np.array([a,b]) for a,b in pairs]
        arches=[p for a,b in pairs[:4] for p in list(candidate_links(a,b))[1:]]
        candidates.extend(sorted(arches,key=lambda p:np.linalg.norm(np.diff(p,axis=0),axis=1).sum()))
        for path in candidates:
                attempts+=1
                if attempts%100==0:print(args.module,'tested local ribs',attempts,flush=True)
                raw=D.outer_envelope(path,1.)
                # A fitted TPU rib may conform to a bone or go around a cap;
                # check the actual trimmed solid, never exempt those pairs.
                rib=raw
                for mask in masks:
                    if bbox_overlap(rib,mask):rib=rib-mask
                rib=rib-lumens
                if rib.volume()<=1e-6:continue
                # Most proposals lose one attachment when fitted around an
                # obstacle. Reject them before an expensive full-mesh union.
                if I.intersection_volume(rib,fragment)<=1e-5:continue
                if I.intersection_volume(rib,main_material)<=1e-5:continue
                hits={k:I.intersection_volume(rib,v) for k,v in obstacles.items()}
                if any(v>.001 for v in hits.values()):continue
                joined=solid+rib
                after=sum(p.volume()>1e-9 for p in joined.decompose())
                if after>=len(now):continue
                found=(joined,dict(path_source_mm=path.tolist(),radius_mm=1.,
                    length_mm=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum()),
                    new_material_mm3=float((joined-solid).volume()),
                    raw_rib_volume_mm3=raw.volume(),shaped_rib_volume_mm3=rib.volume(),
                    obstacle_intersections_mm3=hits,before_components=len(now),after_components=after))
                break
        if found:
            solid,row=found;links.append(row)
            print(args.module,'joined fragment',fragment.volume(),'components',row['before_components'],'->',row['after_components'],flush=True)
        else:
            unresolved.append(dict(volume_mm3=fragment.volume(),bounds_source_mm=list(fragment.bounding_box())))
            print(args.module,'no clear local rib for fragment',fragment.volume(),flush=True)
    out=D.DEST/('unified_'+args.module+'_connected');out.mkdir(exist_ok=True)
    restored,part=local_export(solid,out/'TPU_muscle.stl',origin)
    parts={**meta['parts'],'TPU':part}
    report=dict(meta,scope='static TPU transition-rib candidate; not validated compliance or stage acceptance',
        parts=parts,parent_report=str(parent_file.relative_to(D.frames.PROJECT_ROOT)),
        parent_report_sha256=hashlib.sha256(parent_file.read_bytes()).hexdigest(),
        original_positive_components=original_components,transition_ribs=links,
        tested_candidates=attempts,unresolved_fragments=unresolved,
        added_rib_strength_and_deformation_verified=False,
        repaired_TPU_bone_and_neighbour_intersections_mm3={k:I.intersection_volume(restored,v) for k,v in obstacles.items()},
        stage_pass=False)
    # Parent numerical checks describe a different TPU mesh. Recompute the
    # physical steel and cap checks and explicitly invalidate old wall data.
    report['parent_free_TPU_wall_guard']=report.pop('free_TPU_wall_guard_outside_assigned_cap_housings')
    report['parent_nominal_wall_guard']=report.pop('missing_0_75mm_wall_guard_mm3')
    report['wall_guard_after_ribs_verified']=False
    allparts={key:load_row(row) for key,row in parts.items()}
    report['wire_material_intersections_mm3']={r['head']:{k:I.intersection_volume(
        D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3),s) for k,s in allparts.items()} for r in routes}
    report['material_pair_intersections_mm3']=[dict(a=a,b=b,intersection_mm3=I.intersection_volume(allparts[a],allparts[b]))
        for a,b in itertools.combinations(allparts,2)]
    report['all_part_bone_intersections_mm3']={k:{b:I.intersection_volume(s,obstacles[b]) for b in D.SP.BONES} for k,s in allparts.items()}
    report['bone_intersections_mm3']=report['all_part_bone_intersections_mm3']['TPU']
    report['parent_neighbour_intersections_mm3']=report.pop('unchanged_neighbour_intersections_mm3')
    report['nominal_hand_pull_endpoint_intersections_mm3']={r['head']:{k:I.intersection_volume(
        I._manifold().Manifold.sphere(.3,32).translate(r['ports']['origin']['hand_pull_lead_source_mm']),s)
        for k,s in allparts.items()} for r in routes}
    (out/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(args.module,'mesh',part['mesh_checks_pass'],'ribs',len(links),'unresolved',len(unresolved),flush=True)


if __name__=='__main__':main()
