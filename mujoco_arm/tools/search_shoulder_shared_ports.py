"""Finite search over three shared insertion ports, with bone envelopes.

Preserves the bone seat, cap footprint and three existing active core paths.
Replaces the old common-centre short guide paths, which were not separated.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import build_shoulder_tendon_segments as T
import mesh_geom
import shoulder_pads
import interference as I


def main():
    source=D.DEST/'routes.json';report=json.loads(source.read_text(encoding='utf8'))
    heads=['FJ1468','FJ1513','FJ1467'];by={r['head']:r for r in report['routes']}
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    seat=next(s for s in seats if len(s['heads'])==3)
    c=np.array(seat['centre_source_mm']);n=D.unit(seat['normal_source']);basis=shoulder_pads.pad_basis(n)
    height=seat['lift_mm']+seat['pad_thickness_mm']+3
    bones={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()};tris=[]
    for v in D.SP.BONES.values():
        V,F=mesh_geom.load_obj(str(T.BONE_DIR/v));W,G=mesh_geom.weld(V,F);tris.append(W[G])
    trials=[];selected=None
    for phase in [30,25,35,20,40]:
        offsets=np.array([basis[:,0]*4.5*np.cos(a)+basis[:,1]*4.5*np.sin(a)
                          for a in np.radians(phase+np.array([0,120,240]))])
        # Two existing M3 hole centres; r1.4 crimp seats need >=0.5mm wall.
        bolts=[basis[:,0]*s*7 for s in [-1,1]]
        bolt_wall=min(np.linalg.norm(p-b)-1.4-1.7 for p in offsets for b in bolts)
        if bolt_wall<.5:continue
        for perm in itertools.permutations(range(3)):
            choices={};all_clear=True
            for h,i in zip(heads,perm):
                inner=c+offsets[i]+n*(height+.2);outer=inner+n*seat['cap_thickness_mm'];bend=inner-n*2
                endpoint=np.array(by[h]['core_path_source_mm'][-1])
                path=np.vstack([inner,T.relax_off_bone([bend,endpoint],tris,3.5,step_mm=1,iterations=60,lam=.2)])
                env=D.outer_envelope(path,r=3.2)
                hits={k:I.intersection_volume(env,b) for k,b in bones.items()}
                all_clear=all_clear and all(v<=.001 for v in hits.values())
                full=np.vstack([by[h]['ports']['origin']['hand_pull_lead_source_mm'],
                    by[h]['ports']['origin']['outer_face_source_mm'],
                    by[h]['connectors']['origin']['path_source_mm'],
                    np.array(by[h]['core_path_source_mm'])[1:],path[-2::-1],outer])
                choices[h]=dict(path=path,full=full,inner=inner,outer=outer,offset=offsets[i],hits=hits)
            separations=[dict(a=a,b=b,distance_mm=D.path_distance(choices[a]['full'],choices[b]['full']))
                         for a,b in itertools.combinations(heads,2)]
            distance=min(r['distance_mm'] for r in separations)
            trial=dict(phase_deg=phase,assignment=list(perm),minimum_centreline_distance_mm=distance,
                       bone_envelopes_clear=all_clear,min_bolt_to_crimp_wall_mm=float(bolt_wall))
            trials.append(trial);print('Shared port trial',trial,flush=True)
            if all_clear and distance>=2.45:
                selected=dict(trial=trial,pair_separations=separations);break
        if selected:break
    if selected:
        for h in heads:
            choice=choices[h];row=by[h]
            row['ports']['distal'].update(offset_from_seat_mm=choice['offset'].tolist(),
                inner_face_source_mm=choice['inner'].tolist(),outer_face_source_mm=choice['outer'].tolist())
            row['connectors']['distal'].update(path_source_mm=choice['path'].tolist(),
                outer_envelope_bone_intersection_mm3=choice['hits'],bone_clear_at_nominal_geometry=True)
            row['full_wire_path_source_mm']=choice['full'].tolist()
        report['deltoid_wire_pair_minimum_mm']=selected['pair_separations']
        report['shared_port_assignment_search']=selected
        report['scope']='static unified routes with separated shared insertion ports; NOT manufactured or moving acceptance'
        (D.DEST/'routes_separated.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    (D.DEST/'shared_port_search.json').write_text(json.dumps(dict(trials=trials,selected=selected),indent=2),encoding='utf8')
    print('Separated shared-port candidate found:',bool(selected),flush=True)


if __name__=='__main__':main()
