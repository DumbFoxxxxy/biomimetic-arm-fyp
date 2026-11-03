"""Integral keyed GH nut pockets with removable transverse M2 keepers.

Parallel pocket cheeks react nominal M3 nut rotation; an M2 cross-screw
closes the side entrance for nut and washer. Thread/preload, material and
physical holding strength remain unverified. No muscle or joint relocation.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_bench_mount as B
from open_shoulder_GH_assembly_entries import DEST as PARENT
from build_shoulder_shaft_collar import rod
from build_shoulder_unified_muscle import local_export
from check_shoulder_front_liner_keepers import is_compliant

DEST=A.DESIGN/'s1_GH_assembly_review/nut_keeper_candidate'


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    source=PARENT/'assembly_check.json';parent=json.loads(source.read_text(encoding='utf8'))
    sources=dict(parent['sources']);solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in sources.items()}
    reference=A.G1/'assembly.json';joint=json.loads(reference.read_text(encoding='utf8'))
    G=np.array(joint['centre_source_mm']);n=np.array(joint['neck_axis']);mf=A.I._manifold()
    old=solids['bone:scapula'];body=old;definitions=[]
    for i,entry in enumerate(parent['GH_assembly_entries']['definitions'][:2]):
        u=np.array(entry['direction_source']);v=np.cross(n,u);F=np.column_stack([u,v,n])
        c=np.array(joint['bolt_sites_display_mm'][i])+A.frames.SHIFT
        boss=B.box([-4.6,-5.5,-8.5],[7.2,5.5,-3.0])
        # Different heights accommodate the wider washer above the keyed nut.
        boss-=B.box([-3.35,-2.85,-7.3],[8.,2.85,-4.7])
        boss-=B.box([-3.65,-3.6,-4.8],[8.,3.6,-3.95])
        boss-=rod([0,0,-10],[0,0,1],1.7)
        boss-=rod([4.7,-6.,-5.2],[4.7,6.,-5.2],1.1)
        boss=boss.transform(np.c_[F,c])-mf.Manifold.sphere(joint['socket_inner_radius_mm'],40).translate(G)
        body+=boss
        axis=v;point=c+F@np.array([4.7,-5.9,-5.2])
        shank=rod(point,point+axis*14.,1.)
        head=rod(point-axis*2.,point,1.9)
        screw=shank+head
        near=rod(point,point+axis*.4,2.25)-rod(point-axis*.1,point+axis*.5,1.1)
        far=rod(point+axis*11.4,point+axis*11.8,2.25)-rod(point+axis*11.3,point+axis*11.9,1.1)
        w=A.D.unit(np.cross(np.eye(3)[np.argmin(abs(axis))],axis));H=np.column_stack([w,np.cross(axis,w),axis])
        nut=mf.Manifold.cylinder(1.6,4./np.sqrt(3),4./np.sqrt(3),6)-mf.Manifold.cylinder(1.8,1.1,1.1,48).translate([0,0,-.1])
        nut=nut.transform(np.c_[H,point+axis*11.8])
        hardware={}
        for name,part in dict(screw=screw,near_washer=near,far_washer=far,nut=nut).items():
            key=f'GH_keeper:{i}/{name}'
            actual,meta=local_export(part,DEST/f'keeper_{i}_{name}.stl',c)
            if not meta['mesh_checks_pass']:raise ValueError('Keeper mesh fails '+key)
            solids[key]=actual;sources[key]=meta;hardware[name]=key
        definitions.append(dict(index=i,source_origin_mm=c.tolist(),basis_columns_source=F.tolist(),
                                entry_direction_source=u.tolist(),keyed_nut_channel_width_mm=5.7,
                                nominal_M3_nut_across_flats_mm=5.5,washer_channel_width_mm=7.2,
                                boss_local_bounds_mm=[-4.6,-5.5,-8.5,7.2,5.5,-3.],
                                cross_pin_local_u_n_mm=[4.7,-5.2],cross_bore_diameter_mm=2.2,
                                keeper_axis_source=axis.tolist(),nominal_M2_screw_length_mm=14.,
                                keeper_hardware=hardware,old_bone_overlap_mm3=A.I.intersection_volume(boss,old)))
    if sum(p.volume()>1e-6 for p in body.decompose())!=1:raise ValueError('Key pockets disconnected from engineered scapula')
    origin=np.array(parent['sources']['bone:scapula']['translation_source_mm'])
    actual,meta=local_export(body,DEST/'scapula_keyed_GH_pockets.stl',origin)
    if not meta['mesh_checks_pass']:raise ValueError('Keyed scapula mesh fails')
    sources['bone:scapula']=meta;solids['bone:scapula']=actual
    changes=dict(added_mm3=float((actual-old).volume()),removed_by_export_mm3=float((old-actual).volume()))
    overlaps=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=volume,compliant=is_compliant(a) or is_compliant(b)))
    routes=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';wires=[]
    for r in json.loads(routes.read_text(encoding='utf8'))['routes']:
        wire=A.D.outer_envelope(np.array(r['full_wire_path_source_mm']),.3)
        hits={k:A.I.intersection_volume(wire,s) for k,s in solids.items()}
        wires.append(dict(head=r['head'],hits_over_threshold=[dict(part=k,volume_mm3=z) for k,z in hits.items() if z>.001]))
    data=dict(scope=__doc__,n_instances=len(solids),sources=sources,parent_report=source.relative_to(A.ROOT).as_posix(),
              parent_sha256=A.hashed(source),script_sha256=A.hashed(Path(__file__).resolve()),
              joint_reference_sha256=A.hashed(reference),routes_sha256=A.hashed(routes),
              GH_nut_keepers=dict(definitions=definitions,scapula_changes=changes),
              GH_assembly_entries=parent['GH_assembly_entries'],GH_retention_revision=parent['GH_retention_revision'],GH_lip_relief=parent['GH_lip_relief'],
              free_terminal_constraints_required=parent['free_terminal_constraints_required'],
              intended_body_owners={**parent['intended_body_owners'],**{k:'fixed_context' for k in sources if k.startswith('GH_keeper:')}},
              pair_count_checked=len(solids)*(len(solids)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
              stage_pass=False,actual_drive_range_deg=None,
              not_verified=['pocket key reaction, finite capture and staged installation need independent readback checks',
                            'M2 keeper assembly and real tool/head access, nut reaction and preload',
                            'printed tolerances, material, pocket wall and keeper strength',
                            'full shoulder muscle placement, deformation and physical actuation'])
    (DEST/'assembly_check.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    print('Scapula change',changes,flush=True)
    print('Rigid overlaps',[r for r in overlaps if r['volume_mm3']>.001 and not r['compliant']],flush=True)
    print('New soft overlaps',[r for r in overlaps if r['volume_mm3']>.001 and r['compliant'] and (r['a']=='bone:scapula' or r['b']=='bone:scapula' or r['a'].startswith('GH_keeper:') or r['b'].startswith('GH_keeper:'))],flush=True)
    print('Wire failures',[r for r in wires if r['hits_over_threshold']],flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
