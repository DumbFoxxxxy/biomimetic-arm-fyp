"""Independent split-cap topology, lateral placement and wire/stop review.

Only cap replacement stage is reviewed. Earlier muscle placement is still
blocked by collar parts and must not be declared solved by this local result.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component
import mesh_geom

DEST=A.DESIGN/'s1_teres_split_cap'


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/'tools/build_s1_teres_split_cap.py')==r['script_sha256']
    for key,row in parent['sources'].items():
        assert sha(ROOT/row['file'])==row['sha256']
        if key not in r['retired_keys']:assert row==r['sources'][key]
    solids={key:load_component(row) for key,row in r['sources'].items()}
    topology=[]
    for key in r['new_keys']:
        v,f=mesh_geom.load_obj(str(ROOT/r['sources'][key]['file']))
        _,inv=np.unique(v,axis=0,return_inverse=True);faces=inv[f]
        _,bnd,nm=mesh_geom.edge_census(faces);xyz=v[f]
        deg=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        area0=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        count=sum(piece.volume()>1e-8 for piece in solids[key].decompose())
        assert bnd==nm==deg==area0==0 and count==1,key
        topology.append(dict(key=key,boundary=bnd,nonmanifold=nm,repeated=deg,zero_area=area0,positive_components_over_1e_8_mm3=count))
    later={key for key in solids if key.startswith(('distal:FJ1507/','stop:FJ1507/','wire_active:FJ1507/'))}
    fixed=set(solids)-set(r['new_keys'])-later
    assert len(fixed)==199
    start_keys=sorted(fixed);F=np.array(r['local_basis_columns']);normal=F[:,0];n=F[:,2]
    stages=[]
    for suffix,sign in [('A',-1),('B',1)]:
        key='cap:FJ1507/distal_'+suffix;hits=[]
        for way,values in [('in',np.arange(20.,-.001,-.5)),('out',np.arange(0.,20.001,.5))]:
            for distance in values:
                posed=solids[key].translate(normal*sign*distance)
                for target in fixed:
                    volume=float(A.I.intersection_volume(posed,solids[target]))
                    if volume>.001:hits.append(dict(way=way,offset_mm=float(distance),target=target,volume_mm3=volume))
        stages.append(dict(key=key,sign=sign,direction_source=(normal*sign).tolist(),step_mm=.5,poses=82,
                           fixed_keys=sorted(fixed),hits=hits))
        fixed.add(key)
        print(suffix,'hits',len(hits),flush=True)
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    assert sha(route_file)==r['routes_sha256']
    corridors=[]
    for row in read(route_file)['routes']:
        wire=A.D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
        volumes={key:float(A.I.intersection_volume(wire,solids[key])) for key in r['new_keys']}
        corridors.append(dict(head=row['head'],cap_overlap_mm3=volumes))
    combined=solids[r['new_keys'][0]]+solids[r['new_keys'][1]]
    stop=solids['stop:FJ1507/hex_body']
    bearing=float(A.I.intersection_volume(stop.translate(-n*.08),combined))
    proof=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),topology=topology,
        initial_199_fixed_keys=start_keys,later_installation_keys=sorted(later),stages=stages,
        all_other_206_sources_unchanged=True,six_active_wire_corridors=corridors,
        stop_bearing_probe_shift_mm=.08,stop_bearing_probe_mm3=bearing,
        split_cap_geometry_pass=(r['static_geometry_pass'] and not any(row['hits'] for row in stages) and bearing>.001
             and all(value<=.001 for row in corridors for value in row['cap_overlap_mm3'].values())),
        real_pin_fit_verified=False,strength_verified=False,complete_module_installation_verified=False,
        driven_motion_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',proof)
    print('Split cap pass',proof['split_cap_geometry_pass'],'bearing',bearing,flush=True)
    if not proof['split_cap_geometry_pass']:raise ValueError('Split cap approach needs revision; saved conflicts retained')


if __name__=='__main__':main()
