"""Exact solid counterexample for direct fixed-sheave -> receiver wiring.

Each span departs a genuine root-circle point, never a centre/control point.
Receiver axes are conditionally aligned to that point. An out-of-plane corner
is NOT accepted as a valid smooth wire path even if the finite solid is clear.
"""
from screen_forearm_receivers import receiver_pose
from build_forearm_receivers import *
MOUNT=DEST/'rephased_anchors'

def run():
    screen=json.loads((MOUNT/'mounted_direction_screen.json').read_text());definition=json.loads((BASE/'design.json').read_text());axis=np.array(definition['joint_axis']);E=np.array(definition['joint_centre_display_mm']);down=np.array(definition['zero_forearm_direction']);F0=np.stack([down,np.cross(-axis,down),-axis],axis=1);selection=json.loads((OUT/'elbow_guides/paired_placement/selection.json').read_text())['placements'];reports=[]
    for fid,label in [('FJ1486','flexor'),('FJ1477','extensor')]:
        s=selection[fid];Q=np.array(s['centre_display_mm']);F=np.array(s['frame_columns']);P0=np.array(screen['candidate_eyelets_display_mm'][label]);d=json.loads((DEST/fid/'design.json').read_text());parts={k:load(OUT/'elbow_guides/paired_placement'/(fid+'_'+k+'.obj'))[0] for k in json.loads((OUT/'elbow_guides/template/design.json').read_text())['parts']};rows=[];radius=8.005
        for angle in range(151):
            R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();P=E+R@(P0-E);frame=R@F0;target=Q.copy();converged=False
            for iteration in range(50):
                _,yaw,pitch,D=receiver_pose({},d['standoff_mm'],frame,P,target);p=(F.T@(D-Q))[:2];distance=np.linalg.norm(p);assert distance>radius;perp=np.array([-p[1],p[0]]);t1=radius**2/distance**2*p+radius*np.sqrt(distance**2-radius**2)/distance**2*perp;t2=radius**2/distance**2*p-radius*np.sqrt(distance**2-radius**2)/distance**2*perp
                t=max([t1,t2],key=lambda t:float((p-t)@np.array([t[1],-t[0]])));next_target=Q+F@np.r_[t,0]
                if np.linalg.norm(next_target-target)<1e-7:converged=True;target=next_target;break
                target=next_target
            _,yaw,pitch,D=receiver_pose({},d['standoff_mm'],frame,P,target);wire=rod(target,D,.3,16);hits=[]
            for name,body in parts.items():
                volume=intersection(wire,body)
                if volume>.001:hits.append({'part':name,'volume_mm3':volume})
            out_of_plane=float(np.degrees(np.arcsin(np.clip(abs(unit(D-target)@F[:,2]),0,1))))
            rows.append({'elbow_angle_deg':angle,'actual_root_departure_display_mm':target.tolist(),'conditional_receiver_outlet_mm':D.tolist(),'receiver_pitch_deg':pitch,'out_of_plane_departure_deg':out_of_plane,'fixed_point_converged':converged,'guide_solid_hits':hits})
        reports.append({'muscle':fid,'samples':rows,'guide_collision_pose_count':sum(bool(r['guide_solid_hits']) for r in rows),'minimum_departure_out_of_plane_deg':min(r['out_of_plane_departure_deg'] for r in rows),'maximum_departure_out_of_plane_deg':max(r['out_of_plane_departure_deg'] for r in rows),'all_tangent_point_iterations_converged':all(r['fixed_point_converged'] for r in rows)})
    inputs=[MOUNT/'mounted_direction_screen.json',BASE/'design.json',OUT/'elbow_guides/paired_placement/selection.json']+[DEST/fid/'design.json' for fid in ['FJ1486','FJ1477']]+list((OUT/'elbow_guides/paired_placement').glob('*.obj'))
    result={'scope':'conditional straight output spans from actual fixed-root-circle tangencies, checked against real7-part sheave; not approved steel routing or installed-gimbal path','modules':reports,'angular_step_elbow_deg':1,'C1_steel_tangent_continuity_verified':False,'wire_minimum_bend_and_friction_verified':False,'full_transmission_exists':False,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'frozen':False};(MOUNT/'direct_departure_probe.json').write_text(json.dumps(result,indent=2),encoding='utf8');print([(r['muscle'],r['guide_collision_pose_count'],r['minimum_departure_out_of_plane_deg'],r['maximum_departure_out_of_plane_deg']) for r in reports],flush=True)

if __name__=='__main__':run()
