"""Nominal tool approach review for the late-installed collar group 0.

Straight envelopes stop before head/nut engagement. They do not establish
real socket fit, hand access, tightening torque or grip.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,read,sha,save,load_component,rod,export

DEST=A.DESIGN/'s1_teres_assembly_order'


def main():
    file=DEST/'assembly_check.json';r=read(file);proof=read(DEST/'installation_checks.json')
    assert proof['assembly_sha256']==sha(file) and proof['nominal_sequence_and_capture_geometry_pass']
    solids={k:load_component(row) for k,row in r['sources'].items()}
    frame_file=A.DESIGN/'s1_retained_joint_comparison/shaft_collar/station_60/interface_check.json'
    frame=read(frame_file)['frame_source'];F=np.array(frame['basis_columns']);I=np.array(frame['origin_mm']);axis=F[:,1]
    local=lambda solid:solid.transform(np.c_[F.T,-F.T@I])
    b=np.array(local(solids['collar:bolt_0']).bounding_box()).reshape(2,3)
    nut=np.array(local(solids['collar:nut_0']).bounding_box()).reshape(2,3)
    mid=(b[0]+b[1])/2
    front=I+F@np.array([mid[0],b[1,1]+.1,mid[2]])
    rear=I+F@np.array([mid[0],nut[0,1]-.1,mid[2]])
    driver=rod(front,front+axis*25.,1.5)
    # Compact socket cup with a narrow shank, instead of a full-length wide
    # tube that strikes the cap. This is a custom nominal tool envelope.
    cup=rod(rear-axis*5.,rear,3.6)-rod(rear-axis*5.1,rear+axis,1.7)
    mf=A.I._manifold();hex_frame=np.stack([F[:,0],-F[:,2],axis],axis=1)
    hex_cavity=mf.Manifold.cylinder(3.7,5.7/np.sqrt(3),5.7/np.sqrt(3),6).transform(
        np.c_[hex_frame,rear-axis*3.2])
    socket=(cup-hex_cavity)+rod(rear-axis*5.,rear-axis*4.2,3.6)+rod(rear-axis*20.,rear-axis*4.5,1.5)
    result=[];folder=DEST/'tool_reference';folder.mkdir(exist_ok=True);sources={}
    for key,tool in [('tool:collar_front_driver',driver),('tool:collar_rear_socket',socket)]:
        hits=[]
        for target,solid in solids.items():
            volume=float(A.I.intersection_volume(tool,solid))
            if volume>.001:hits.append(dict(target=target,volume_mm3=volume))
        _,sources[key]=export(tool,folder/(key.split(':')[1]+'.obj'),I)
        result.append(dict(tool=key,hits=hits))
    approach=[]
    for way,distances in [('in',np.arange(4.,-.001,-.5)),('out',np.arange(0.,4.001,.5))]:
        for distance in distances:
            for target,solid in solids.items():
                volume=float(A.I.intersection_volume(socket.translate(-axis*distance),solid))
                if volume>.001:approach.append(dict(way=way,offset_mm=float(distance),target=target,volume_mm3=volume))
    save(DEST/'collar_tool_checks.json',dict(scope=__doc__,assembly_sha256=sha(file),script_sha256=sha(Path(__file__)),
        frame_sha256=sha(frame_file),tool_sources=sources,tool_checks=result,
        socket_dimensions=dict(cup_OD_mm=7.2,cup_length_mm=5.,shank_OD_mm=3.,nominal_hex_AF_mm=5.7),
        local_socket_entry_mm=4.,entry_step_mm=.5,entry_poses=18,entry_hits=approach,
        simultaneous_tools_overlap_mm3=float(A.I.intersection_volume(driver,socket)),
        nominal_straight_tool_approach_pass=not any(row['hits'] for row in result) and not approach,
        actual_tool_engagement_verified=False,operator_hand_space_verified=False,stage_pass=False))
    print(result,'Entry hits',approach[:4],len(approach),flush=True)


if __name__=='__main__':main()
