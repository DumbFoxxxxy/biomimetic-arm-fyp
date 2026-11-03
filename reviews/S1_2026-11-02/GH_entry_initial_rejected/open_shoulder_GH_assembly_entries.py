"""Local staged assembly entrances in the engineered scapula GH connection.

Keep anatomical shaft/blade outside GH neighbourhood, humerus and muscles.
Continuous swept keepouts open existing nut/washer pockets and cap/bolt access.
The resulting pocket walls and retained shell require strength/print testing.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_bench_mount as B
from relieve_shoulder_GH_retaining_caps import DEST as PARENT
from build_shoulder_unified_muscle import local_export
from check_shoulder_front_liner_keepers import is_compliant

DEST = A.DESIGN/'s1_GH_assembly_review/assembly_entries'


def swept(body, direction, length, allowance=.1):
    helper = np.eye(3)[np.argmin(abs(direction))]
    x = A.D.unit(np.cross(helper, direction))
    F = np.column_stack([x, np.cross(direction, x), direction])
    prism = B.box([-allowance]*3, [allowance,allowance,length+allowance])
    return body.minkowski_sum(prism.transform(np.c_[F,np.zeros(3)]))


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    source=PARENT/'assembly_check.json';parent=json.loads(source.read_text(encoding='utf8'))
    sources=dict(parent['sources'])
    solids={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in sources.items()}
    jointfile=A.G1/'assembly.json';joint=json.loads(jointfile.read_text(encoding='utf8'))
    n=np.array(joint['neck_axis']);G=np.array(joint['centre_source_mm'])
    x=A.D.unit(np.cross(n,[0.,-1.,0.]));y=np.cross(n,x)
    old=solids['bone:scapula']; revised=old; definitions=[]
    for i,angle in [(0,240.),(1,60.)]:
        d=x*np.cos(np.radians(angle))+y*np.sin(np.radians(angle))
        keys=[f'GH_nut_{i}',f'GH_bottom_washer_{i}']
        body=solids[keys[0]]+solids[keys[1]]
        keepout=swept(body,d,25.)
        removal=revised^keepout; revised-=keepout
        definitions.append(dict(name=f'nut_and_bottom_washer_{i}',moving_keys=keys,direction_source=d.tolist(),
                                maximum_translation_mm=25.,removed_mm3=float(removal.volume())))
    for label,keys,length in [('cap_and_top_washer_1',['GH_cap_1','GH_top_washer_1'],25.),
                               ('bolt_1',['GH_bolt_1'],15.)]:
        body=A.I._manifold().Manifold.batch_boolean([solids[k] for k in keys],A.I._manifold().OpType.Add)
        keepout=swept(body,n,length); removal=revised^keepout;revised-=keepout
        definitions.append(dict(name=label,moving_keys=keys,direction_source=n.tolist(),
                                maximum_translation_mm=length,removed_mm3=float(removal.volume())))
    if sum(p.volume()>1e-6 for p in revised.decompose())!=1:
        raise ValueError('GH entry opening disconnects engineered scapula')
    actual,meta=local_export(revised.as_original().simplify(.001),DEST/'scapula_GH_entries.stl',G)
    if not meta['mesh_checks_pass']:raise ValueError('Entry scapula mesh readback fails')
    removal=old-actual;outside=removal-A.I._manifold().Manifold.sphere(55,48).translate(G)
    if outside.volume()>.001:raise ValueError('Entry modification extends outside GH engineering neighbourhood')
    sources['bone:scapula']=meta;solids['bone:scapula']=actual
    overlaps=[]
    for a,b in itertools.combinations(solids,2):
        volume=A.I.intersection_volume(solids[a],solids[b])
        if volume>1e-9:overlaps.append(dict(a=a,b=b,volume_mm3=volume,compliant=is_compliant(a) or is_compliant(b)))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';wires=[]
    for row in json.loads(routefile.read_text(encoding='utf8'))['routes']:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        hits={k:A.I.intersection_volume(wire,s) for k,s in solids.items()}
        wires.append(dict(head=row['head'],volumes_mm3=hits,hits_over_threshold=[dict(part=k,volume_mm3=v) for k,v in hits.items() if v>.001]))
    data=dict(scope=__doc__,n_instances=len(solids),sources=sources,
              parent_report=source.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(source),
              script_sha256=A.hashed(Path(__file__).resolve()),joint_reference_sha256=A.hashed(jointfile),
              routes_sha256=A.hashed(routefile),
              generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                  [A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',A.ROOT/'mujoco_arm/tools/build_shoulder_bench_mount.py']},
              GH_assembly_entries=dict(continuous_sweep_allowance_mm=.1,definitions=definitions,
                  changed_key='bone:scapula',removed_mm3=float(removal.volume()),added_by_export_mm3=float((actual-old).volume()),
                  removed_outside_source_GH_radius55mm_mm3=float(outside.volume()),positive_components=1),
              GH_retention_revision=parent['GH_retention_revision'],GH_lip_relief=parent['GH_lip_relief'],
              pair_count_checked=len(solids)*(len(solids)-1)//2,material_intersections=overlaps,actuation_wire_checks=wires,
              free_terminal_constraints_required=parent['free_terminal_constraints_required'],
              intended_body_owners=parent['intended_body_owners'],
              scapula_insertion_relief=parent['scapula_insertion_relief'],keeper_definitions=parent['keeper_definitions'],
              stage_pass=False,actual_drive_range_deg=None,
              not_verified=['dense staged assembly and capture must be rechecked after export',
                            'nut antirotation and reaction faces, tool envelopes and temporary holding',
                            'remaining pocket/retained shell strength, printed tolerances and preload',
                            'full arm motion, muscle deformation and physical drive'])
    (DEST/'assembly_check.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    print('Scapula entry modification',data['GH_assembly_entries'],flush=True)
    print('Rigid overlaps',[r for r in overlaps if r['volume_mm3']>.001 and not r['compliant']],flush=True)
    print('Wire failures',sum(bool(r['hits_over_threshold']) for r in wires),flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    run()
