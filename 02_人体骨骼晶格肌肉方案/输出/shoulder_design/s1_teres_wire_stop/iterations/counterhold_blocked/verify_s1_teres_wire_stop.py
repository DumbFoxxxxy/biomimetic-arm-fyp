"""Independent teres stop bearing, wire entry, placement and paired-tool check.

First-touch is a geometry state, never a measured grip or torque capacity.
Counterhold and driver shapes are nominal approach envelopes, not real tools.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,sha,read,save,load_component
import mesh_geom

DEST=A.DESIGN/'s1_teres_wire_stop'


def overlap(a,b):return float(A.I.intersection_volume(a,b))


def main():
    file=DEST/'assembly_check.json';r=read(file);parent=read(ROOT/r['parent_report'])
    assert sha(ROOT/r['parent_report'])==r['parent_sha256']
    assert sha(ROOT/r['parent_report'].replace('assembly_check.json','verification_checks.json'))==r['parent_verification_sha256']
    assert sha(ROOT/'tools/build_s1_teres_wire_stop.py')==r['script_sha256']
    assert len(parent['sources'])==204 and len(r['sources'])==207
    for key,row in parent['sources'].items():assert r['sources'][key]==row and sha(ROOT/row['file'])==row['sha256']
    solids={key:load_component(row) for key,row in r['sources'].items()}
    topology=[]
    for key in r['new_keys']:
        path=ROOT/r['sources'][key]['file'];v,f=mesh_geom.load_obj(str(path))
        _,index=np.unique(v,axis=0,return_inverse=True);faces=index[f]
        _,boundary,nonmanifold=mesh_geom.edge_census(faces)
        repeated=int(np.any(np.diff(np.sort(faces,axis=1),axis=1)==0,axis=1).sum())
        xyz=v[f];zero_area=int((np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)==0).sum())
        positive=sum(piece.volume()>1e-8 for piece in solids[key].decompose())
        topology.append(dict(key=key,boundary=boundary,nonmanifold=nonmanifold,repeated=repeated,zero_area=zero_area,
                             positive_components_over_1e_8_mm3=positive))
        assert boundary==nonmanifold==repeated==zero_area==0 and positive==1,key
    I=np.array(r['local_origin_source_mm']);F=np.array(r['local_basis_columns']);x,y,n=F.T
    body_key='stop:FJ1507/hex_body';screw_key='stop:FJ1507/M2x3';wire_key='wire_active:FJ1507/reference'
    body=solids[body_key];screw=solids[screw_key];wire=solids[wire_key];cap=solids['cap:FJ1507/distal']
    bearing=overlap(body.translate(-n*.08),cap)
    screw_contact=overlap(screw.translate(-x*.05),wire)
    assert bearing>.001 and screw_contact>.001
    # The radial screw is backed off 0.3 mm during axial feeding. The full
    # reference wire is retained in every placement sample, not excluded.
    fixed={key:value for key,value in solids.items() if key not in [body_key,screw_key]}
    placement=[]
    offsets=np.arange(12.,-.001,-.5)
    for way,values in [('in',offsets),('out',offsets[::-1])]:
        for distance in values:
            for key,posed in [(body_key,body.translate(n*distance)),
                              (screw_key,screw.translate(n*distance+x*.3))]:
                for target,obstacle in fixed.items():
                    v=overlap(posed,obstacle)
                    if v>.001:placement.append(dict(way=way,offset_mm=float(distance),moving=key,target=target,volume_mm3=v))
    tightening=[]
    for way,values in [('in',[.3,.2,.1,0.]),('out',[0.,.1,.2,.3])]:
        for distance in values:
            for target,obstacle in {**fixed,body_key:body}.items():
                v=overlap(screw.translate(x*distance),obstacle)
                if v>.001:tightening.append(dict(way=way,offset_mm=distance,target=target,volume_mm3=v))
    # A fork approaches from local -X and bears against the local Y flats.
    # Real spanner dimensions and engagement are still unverified.
    mf=A.I._manifold()
    box=lambda lo,hi:mf.Manifold.cube(np.array(hi)-lo).translate(lo)
    fork=(box([-14.,-4.15,5.3],[-6.,4.15,6.8])+
          box([-8.,2.65,5.3],[1.,4.15,6.8])+
          box([-8.,-4.15,5.3],[1.,-2.65,6.8]))
    fork_world=fork.transform(np.c_[F,I])
    driver=rod(I+F@np.array([3.4,0.,6.05]),I+F@np.array([25.,0.,6.05]),.8)
    tool_records=[]
    for label,tool in [('driver',driver),('counterhold_final',fork_world)]:
        hits=[]
        for key,solid in solids.items():
            v=overlap(tool,solid)
            if v>.001:hits.append(dict(key=key,volume_mm3=v))
        tool_records.append(dict(tool=label,hits=hits))
    fork_path=[]
    for way,values in [('in',np.arange(20.,-.001,-.5)),('out',np.arange(0.,20.001,.5))]:
        for distance in values:
            tool=fork_world.translate(-x*distance)
            for key,solid in solids.items():
                v=overlap(tool,solid)
                if v>.001:fork_path.append(dict(way=way,offset_mm=float(distance),key=key,volume_mm3=v))
    dual_tool_overlap=overlap(driver,fork_world)
    jaw_probe=[]
    for sign in [-1,1]:
        jaw=(box([-1.,2.45 if sign>0 else -2.8,5.3],[1.,2.8 if sign>0 else -2.45,6.8])).transform(np.c_[F,I])
        jaw_probe.append(overlap(jaw,body))
    assert min(jaw_probe)>.001
    result=dict(scope=__doc__,assembly_sha256=sha(file),verifier_sha256=sha(Path(__file__)),parent_204_sources_unchanged=True,
        topology=topology,bearing_shift_mm=.08,bearing_probe_mm3=bearing,
        screw_first_contact_shift_mm=.05,screw_wire_contact_probe_mm3=screw_contact,
        stop_to_cap_min_gap_capped_2mm=float(A.I.gap(body,cap,2.)),
        axial_feed=dict(step_mm=.5,poses=50,screw_backoff_mm=.3,hits=placement),
        radial_screw_to_first_touch=dict(step_mm=.1,poses=8,hits=tightening),
        tools=tool_records,counterhold_approach=dict(step_mm=.5,poses=82,hits=fork_path),
        simultaneous_tools_overlap_mm3=dual_tool_overlap,counterhold_flat_probe_mm3=jaw_probe,
        geometry_checks_pass=(r['static_geometry_pass'] and not placement and not tightening
             and not fork_path and not any(row['hits'] for row in tool_records) and dual_tool_overlap<=.001),
        actual_tool_engagement_verified=False,physical_grip_verified=False,grip_force_N=None,tightening_torque_Nm=None,
        complete_module_installation_verified=False,TPU_deformation_verified=False,stage_pass=False,actual_drive_range_deg=None)
    save(DEST/'verification_checks.json',result)
    print(dict(pass_geometry=result['geometry_checks_pass'],bearing=bearing,wire_probe=screw_contact,
        axial_hits=len(placement),radial_hits=len(tightening),tools=tool_records,fork_path_hits=len(fork_path)),flush=True)
    if not result['geometry_checks_pass']:raise ValueError('Stop/tool layout needs revision; see recorded conflicts')


if __name__=='__main__':main()
