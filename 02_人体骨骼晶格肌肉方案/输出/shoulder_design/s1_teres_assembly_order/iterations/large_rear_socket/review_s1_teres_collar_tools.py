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
    socket=rod(rear-axis*20.,rear,4.25)-rod(rear-axis*21.,rear+axis,1.7)
    result=[];folder=DEST/'tool_reference';folder.mkdir(exist_ok=True);sources={}
    for key,tool in [('tool:collar_front_driver',driver),('tool:collar_rear_socket',socket)]:
        hits=[]
        for target,solid in solids.items():
            volume=float(A.I.intersection_volume(tool,solid))
            if volume>.001:hits.append(dict(target=target,volume_mm3=volume))
        _,sources[key]=export(tool,folder/(key.split(':')[1]+'.obj'),I)
        result.append(dict(tool=key,hits=hits))
    save(DEST/'collar_tool_checks.json',dict(scope=__doc__,assembly_sha256=sha(file),script_sha256=sha(Path(__file__)),
        frame_sha256=sha(frame_file),tool_sources=sources,tool_checks=result,
        nominal_straight_tool_approach_pass=not any(row['hits'] for row in result),
        actual_tool_engagement_verified=False,operator_hand_space_verified=False,stage_pass=False))
    print(result,flush=True)


if __name__=='__main__':main()
