"""Two-piece shaft-side muscle anchor comparison on actual engineered humerus.

Includes contoured split TPU liners, two M3x16 screw envelopes, nuts/washers,
and integral tendon eyes. It is an interface candidate, not a connected TPU
muscle design or a friction/strength/drive acceptance result.
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
from build_shoulder_isolated_retained_joint import export_checked
from build_shoulder_tendon_segments import rod
from analyze_shoulder_local_actuation import tension_moments, positive_spanning

DEST=A.DESIGN/'s1_retained_joint_comparison/shaft_collar'
SOURCE=A.DEST/'GH_only_candidate/assembly_check.json'


def box(lo,hi):
    return A.I._manifold().Manifold.cube(np.array(hi)-lo).translate(lo)


def screw(site,axis,length=16.):
    z=A.D.unit(axis);x=A.D.unit(np.cross(z,[0,0,1]));y=np.cross(z,x)
    head=5.25
    profile=np.array([[0,head-length],[1.5,head-length],[1.5,head],
                      [2.75,head],[2.75,head+3],[0,head+3]])
    mf=A.I._manifold()
    return mf.Manifold.revolve(mf.CrossSection([profile]),64).transform(
        np.c_[np.stack([x,y,z],axis=1),site])


def build(station,folder,source):
    folder.mkdir(parents=True,exist_ok=True)
    data=json.loads(source.read_text(encoding='utf8'))
    row=data['sources']['bone:humerus']
    hum,_,_=A.load_checked(A.ROOT/row['file'],row['translation_source_mm'],row['sha256'])
    routes_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(routes_file.read_text(encoding='utf8'))['routes']
    byhead={r['head']:r for r in routes}
    centre=A.frames.joint_centres('source')['GH'];z=A.frames.gh_probe_axes()['longitudinal']
    radial=np.array(byhead['FJ1500']['ports']['distal']['normal_source'])
    y=A.D.unit(radial-z*(radial@z));x=np.cross(y,z);F=np.stack([x,y,z],axis=1)
    local=hum.transform(np.c_[F.T,-F.T@centre])
    half_length=5.;lower=station-half_length;upper=station+half_length
    sections=[local.slice(float(h)) for h in np.arange(lower,upper+.01,.5)]
    mf=A.I._manifold();profile=mf.CrossSection.batch_boolean(sections,mf.OpType.Add).hull()
    inner=profile.offset(1.,circular_segments=48)
    outer=inner.offset(2.5,circular_segments=48)
    b=np.array(outer.bounds());cx=(b[0]+b[2])/2;cy=(b[1]+b[3])/2
    cavity=inner.extrude(upper-lower+4).translate([0,0,lower-2])
    ring=outer.extrude(upper-lower).translate([0,0,lower])-cavity
    liner=inner.extrude(upper-lower-.4).translate([0,0,lower+.2])-local
    # A split liner permits nominal closing at the seam. Its inner surface is
    # exact anatomical contact, not an asserted positive manufacturing gap.
    centre_to_source=lambda s:s.transform(np.c_[F,centre])
    parts={};anchors={};clips={}
    bolt_x=[float(b[0]-2.),float(b[2]+2.)]
    for sign,label,head in [(1,'upper','FJ1500'),(-1,'lower','FJ1504')]:
        clip=box([-100,cy+.5 if sign>0 else -100,lower-10],
                 [100,100 if sign>0 else cy-.5,upper+10])
        clips[label]=clip
        jaw=ring
        for bx in bolt_x:
            jaw=jaw+rod([bx,cy+sign*.5,station],[bx,cy+sign*4.5,station],3.7)
            jaw=jaw+rod([cx,cy+sign*2.5,station],[bx,cy+sign*2.5,station],2.5)
        eyey=float((b[3]+1.8) if sign>0 else (b[1]-1.8))
        eye=np.array([cx,eyey,station])
        jaw=jaw+rod(eye-[1.5,0,0],eye+[1.5,0,0],3.2)
        jaw=jaw+rod([cx,cy,station],eye,2.)
        jaw=(jaw-cavity)^clip
        jaw=jaw-rod(eye-[4,0,0],eye+[4,0,0],1.1)
        for bx in bolt_x:
            jaw=jaw-rod([bx,cy-15,station],[bx,cy+15,station],1.7)
            # The bridge capsules extend beyond the flat ear faces. Make real
            # washer/nut seating pockets instead of excluding their collisions.
            jaw=jaw-rod([bx,cy+sign*4.5,station],[bx,cy+sign*14,station],3.75)
        parts['jaw_'+label]=centre_to_source(jaw)
        parts['liner_'+label]=centre_to_source(liner^clip)
        anchors[head]=dict(eye_centre_source_mm=(F@eye+centre).tolist(),
            bore_axis_source=x.tolist(),eye_bore_radius_mm=1.1,
            proposed_muscle_terminal_source_mm=byhead[head]['ports']['distal']['inner_face_source_mm'],
            straight_link_length_mm=float(np.linalg.norm(F@eye+centre-np.array(byhead[head]['ports']['distal']['inner_face_source_mm']))),
            note='eye only; short tendon, termination and deformable muscle placement NOT designed')
    for index,bx in enumerate(bolt_x):
        site=np.array([bx,cy,station]);axis=np.array([0.,1.,0.]);u=np.array([1.,0.,0.]);v=np.cross(axis,u)
        hardware={
            'bolt':screw(site,axis),
            'top_washer':rod(site+axis*4.65,site+axis*5.15,3.5)-rod(site+axis*4,site+axis*6,1.7),
            'bottom_washer':rod(site-axis*5.15,site-axis*4.65,3.5)-rod(site-axis*6,site-axis*4,1.7),
            'nut':mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(
                np.c_[np.stack([u,v,axis],axis=1),site-axis*7.65])-
                rod(site-axis*9,site-axis*4.5,1.6)}
        for key,solid in hardware.items():parts[f'{key}_{index}']=centre_to_source(solid)
    actual={};meta={}
    for key,solid in parts.items():
        prepared=solid.as_original().simplify(.001)
        actual[key],meta[key]=export_checked(prepared,folder/(key+'.stl'))
        meta[key]['pre_export_simplification_mm']=.001
        meta[key]['pre_export_added_mm3']=float((prepared-solid).volume())
        meta[key]['pre_export_removed_mm3']=float((solid-prepared).volume())
    contacts=[]
    for a,bkey in itertools.combinations(actual,2):
        volume=A.I.intersection_volume(actual[a],actual[bkey])
        if volume>1e-9:contacts.append(dict(a=a,b=bkey,volume_mm3=volume))
    bone_intersections={k:A.I.intersection_volume(v,hum) for k,v in actual.items()}
    # This is only a neutral endpoint force-direction proxy. No unsupported
    # endpoint pinning or rigid muscle deformation is introduced.
    origins=np.array([r['ports']['origin']['inner_face_source_mm'] for r in routes])
    insertions=np.array([anchors[r['head']]['eye_centre_source_mm'] if r['head'] in anchors
                         else r['ports']['distal']['inner_face_source_mm'] for r in routes])
    matrix,lengths=tension_moments(origins,insertions,centre)
    links=[]
    for head,anchor in anchors.items():
        start=np.array(anchor['proposed_muscle_terminal_source_mm']);end=np.array(anchor['eye_centre_source_mm'])
        # Deliberately incomplete reference chord: check material obstacles, do
        # not call it an actual tendon route or attach the original cap to it.
        probe=rod(start,end,.3)
        hits=[dict(part=k,volume_mm3=A.I.intersection_volume(probe,s)) for k,s in actual.items()]
        hits.append(dict(part='bone:humerus',volume_mm3=A.I.intersection_volume(probe,hum)))
        links.append(dict(head=head,intersections=hits))
    result=dict(scope=__doc__,station_from_GH_along_anatomical_axis_mm=station,
        dimensions=dict(axial_length_mm=10.,liner_envelope_offset_mm=1.,rigid_wall_mm=2.5,
                        split_gap_mm=1.,M3_clamp_screw_length_mm=16.,sampling_step_mm=.5),
        frame_source=dict(origin_mm=centre.tolist(),basis_columns=F.tolist()),
        assembly_source=source.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(source),
        routes_sha256=A.hashed(routes_file),script_sha256=A.hashed(Path(__file__).resolve()),
        sources=meta,anchors=anchors,part_intersections_mm3=contacts,
        bone_intersections_mm3=bone_intersections,reference_link_probes=links,
        endpoint_proxy=dict(matrix_mm=matrix.tolist(),chord_lengths_mm=lengths.tolist(),
                            torque_cone=positive_spanning(matrix)),
        not_verified=['clamp closing, frictional torque capacity, liner compression and strength',
            'actual short tendons, termination and new flexible muscle ends',
            'motion and combined paths, fixture access, printed tolerances',
            'endpoint force proxy does not model curved tendons or distributed TPU forces'],
        stage_pass=False,actual_drive_range_deg=None)
    (folder/'interface_check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print('Station',station,'parts',len(meta),'pair overlaps >.001',
          [r for r in contacts if r['volume_mm3']>.001],
          'bone >.001',{k:v for k,v in bone_intersections.items() if v>.001},
          'proxy span',result['endpoint_proxy']['torque_cone']['strict_positive_span'],flush=True)
    return result


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--station',type=float,default=50.)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    build(args.station,args.output_dir or DEST/('station_'+format(args.station,'g')),SOURCE)
