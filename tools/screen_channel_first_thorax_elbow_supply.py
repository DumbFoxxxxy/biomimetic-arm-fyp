"""Neutral rest-shape routes between real ferrule seats; no interior fixed via points."""
from itertools import product
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component,rod
from build_channel_first_inventory import provenance


def curve(controls,count):
    a,b,c,d=np.array(controls);t=np.linspace(0.,1.,count)[:,None]
    path=(1-t)**3*a+3*(1-t)**2*t*b+3*(1-t)*t*t*c+t**3*d
    velocity=3*(1-t)**2*(b-a)+6*(1-t)*t*(c-b)+3*t*t*(d-c)
    acceleration=6*(1-t)*(c-2*b+a)+6*t*(d-2*c+b)
    speed=np.linalg.norm(velocity,axis=1)
    if min(speed)<1e-8:return path,0.
    curvature=np.linalg.norm(np.cross(velocity,acceleration),axis=1)/speed**3
    return path,float(1./max(float(curvature.max()),1e-12))


def flat_envelope(path,radius):
    """Flat physical seat faces; no hemispherical cap intruding below seat."""
    mf=A.I._manifold()
    pieces=[rod(a,b,radius) for a,b in zip(path,path[1:])]
    pieces += [mf.Manifold.sphere(radius,24).translate(p) for p in path[1:-1]]
    return mf.Manifold.batch_boolean(pieces,mf.OpType.Add)


def main():
    base=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first'
    assembly_file=base/'thorax_reaction_mount/assembly.json';assembly=read(assembly_file)
    arm_file=base/'elbow_supply_interface/route_endpoints.json';arm=read(arm_file)
    folder=base/'thorax_elbow_supply';folder.mkdir(exist_ok=True)
    keys=list(assembly['sources']);solids=[load_component(assembly['sources'][k]) for k in keys]
    bounds=np.array([s.bounding_box() for s in solids]);accepted={};routes=[]
    def hits(shape):
        b=np.array(shape.bounding_box());mask=np.all(b[:3]<=bounds[:,3:]+1e-8,axis=1)&np.all(bounds[:,:3]<=b[3:]+1e-8,axis=1)
        result=[]
        for i in np.flatnonzero(mask):
            v=float(A.I.intersection_volume(shape,solids[i]))
            if v>.001:result.append(dict(target=keys[i],overlap_mm3=v))
        return result
    for chest in assembly['ports']:
        head=chest['id'];distal=next(r for r in arm['routes'] if r['id']==head)
        q=np.array(chest['ferrule_outlet_source_mm']);u=np.array(chest['outward_axis_source'])
        end=np.array(distal['physical_ferrule_outlet_source_mm']);v=np.array(distal['outward_axis_source'])
        start_seat=np.array(chest['sheath_reaction_seat_source_mm']);end_seat=np.array(distal['physical_sheath_reaction_seat_source_mm'])
        #8mm straight neutral tails clear rigid retention plates; these outer
        #nodes must be elastic/free later, not extra clamped supports.
        proposals=[]
        for first,last in product((40.,70.,100.,140.,200.),repeat=2):
            controls=np.array([q+u*8.,q+u*(8.+first),end+v*(8.+last),end+v*8.])
            p,r=curve(controls,81);length=float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum())+24.
            if r>=8.:proposals.append(dict(handles_mm=[first,last],controls_source_mm=controls.tolist(),coarse_min_radius_mm=r,length_mm=length))
        proposals.sort(key=lambda r:r['length_mm']);trials=[];selected=None
        for row in proposals:
            p,r=curve(row['controls_source_mm'],81);path=np.vstack([start_seat,q,p,end,end_seat])
            envelope=flat_envelope(path,1.5);collision=hits(envelope)
            trial=dict(**row,coarse_hits=collision);trials.append(trial)
            if collision:continue
            gaps={h:float(A.I.gap(envelope,s,.5)) for h,s in accepted.items()}
            trial['coarse_other_supply_gap_mm_capped_half']=gaps
            if min(gaps.values(),default=.5)<.5:continue
            p,r=curve(row['controls_source_mm'],801);path=np.vstack([start_seat,q,p,end,end_seat])
            envelope=flat_envelope(path,1.5);fine=hits(envelope)
            gaps={h:float(A.I.gap(envelope,s,.5)) for h,s in accepted.items()}
            trial.update(fine_min_radius_mm=r,fine_hits=fine,fine_other_supply_gap_mm_capped_half=gaps)
            if fine or r<8. or min(gaps.values(),default=.5)<.5:continue
            steel=flat_envelope(path,.3);steel_hits=hits(steel)
            trial['steel_hits']=steel_hits
            if steel_hits:continue
            accepted[head]=envelope
            selected=dict(**row,path_source_mm=path.tolist(),fine_min_radius_mm=r,
                sample_count=len(path),maximum_chord_mm=float(np.linalg.norm(np.diff(p,axis=0),axis=1).max()),
                geometry_status='neutral_outer_envelope_and_steel_clear',end_supports_only=[
                    dict(owner='thorax',seat_source_mm=start_seat.tolist(),outward_axis_source=u.tolist()),
                    dict(owner='humerus',seat_source_mm=end_seat.tolist(),outward_axis_source=v.tolist())],
                interior_fixed_guides=[],rest_shape_only=True,
                material_and_support_note='Bezier controls define candidate unloaded curvature only. All intermediate nodes including8mm straight tails must move under sheath bending/compression, cable load and contact. End retention,4mm socket tangent law and external reaction forces unverified.')
            break
        routes.append(dict(id=head,candidate_count=len(proposals),tested=len(trials),trials=trials,selected=selected))
        print(head,'candidates',len(proposals),'tested',len(trials),'neutral candidate',selected is not None,flush=True)
    save(folder/'route_screen.json',dict(assembly=provenance(assembly_file),arm_ports=provenance(arm_file),routes=routes,
        generator=provenance(ROOT/'tools/screen_channel_first_thorax_elbow_supply.py'),
        outer_radius_mm=1.5,wire_radius_mm=.3,minimum_bend_radius_exploration_mm=8.,minimum_pair_gap_mm=.5,
        threshold_mm3=.001,coarse_samples=81,fine_samples=801,rigid_obstacles=len(solids),
        scope='Neutral rest-shape geometry only, not installed full supply, motor or force validation.8mm bend criterion uncalibrated. Flat tube ends meet actual compression seats. No backboard. No excluded rigid pairs, cap lumens treated as real holes. Original14 fixed anatomy overlaps not retested by tube screen. Body motion, external shear/bending/compression, friction and actual guides pending.'))


if __name__=='__main__':main()
