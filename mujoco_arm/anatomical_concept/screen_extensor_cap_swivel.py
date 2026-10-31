"""Actual split-trunnion rotor geometry and conditional cap orientation screen.

No fork or pin mount is invented by this screen. A future real fork and axle
must ground the rotor; this is not an installed mechanism or equilibrium proof.
"""
from build_direct_muscle_yoke import *
from scipy.spatial.transform import Rotation
from build_forearm_receivers import bolt

def run():
    assert EXTENSOR_REANCHOR
    fid='FJ1486' if '--flexor' in sys.argv else 'FJ1477';other_fid='FJ1477' if fid=='FJ1486' else 'FJ1486';rotor_key='flexor_cap_rotor' if fid=='FJ1486' else 'extensor_cap_rotor'
    dest=DEST/('cap_swivel_flexor' if fid=='FJ1486' else 'cap_swivel_study');dest.mkdir(exist_ok=True);a=json.loads((DEST/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};offset=np.array(next(p['translation_display_mm'] for p in a['placements'] if p['id']==fid));d=json.loads((muscle_folder(fid)/'design.json').read_text());F=np.array(d['basis_columns_display']);z=F[:,2];A=np.array(d['pad_centres_display_mm'][0])+offset;B=np.array(d['pad_centres_display_mm'][1])+offset;C=B+z*4.15;j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);axis=np.array(j['joint_axis']);trunnion_axis=unit(axis-z*(axis@z));cap=load(muscle_folder(fid)/'distal_rigid_cap.obj')[0].translate(offset);rotor=cap
    for sign in [-1,1]:
        direction=trunnion_axis*sign
        rotor=rotor+rod(C+direction*10,C+direction*15,1.8,32)
        rotor=rotor-rod(C+direction*11,C+direction*16,1.1,32)
    meta=write(dest/(rotor_key+'.stl'),rotor);assert len([s for s in rotor.decompose() if s.volume()>.001])==1
    (dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[meta]},indent=2),encoding='utf8')
    moving={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if k.startswith(fid+'_distal_') or k==fid+'_active_wire_crimp_barrel'};moving[rotor_key]=rotor
    fixed={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items() if r['kind']!='TPU' and k not in moving and k!='direct_muscle_yoke'}
    fixed[other_fid+'_distal_rigid_cap']=load(muscle_folder(other_fid)/'distal_rigid_cap.obj')[0].translate(next(p['translation_display_mm'] for p in a['placements'] if p['id']==other_fid));rows[other_fid+'_distal_rigid_cap']={'owner':'elbow'}
    sourcefiles=[DEST/'assembly.json',BASE/'design.json',muscle_folder(fid)/'design.json',muscle_folder(fid)/'distal_rigid_cap.obj',muscle_folder(other_fid)/'distal_rigid_cap.obj']+[ROOT/r['file'] for r in a['instances']]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles}
    local_hits=[]
    for key,s in moving.items():
        if key==rotor_key:continue
        v=intersection(rotor,s)
        if v>.001:local_hits.append({'part':key,'volume_mm3':v})
    records=[]
    for angle in range(151):
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();centre=E+R@(C-E);target=unit(centre-A);local_target=R.T@target;u=unit(z-trunnion_axis*(z@trunnion_axis));v=unit(local_target-trunnion_axis*(local_target@trunnion_axis));psi=float(np.arctan2(trunnion_axis@np.cross(u,v),u@v));Q=R@Rotation.from_rotvec(trunnion_axis*psi).as_matrix();T=np.c_[Q,centre-Q@C];hits=[]
        for key,s in moving.items():
            candidate=s.transform(T)
            for other,body in fixed.items():
                target_body=body.transform(np.c_[R,E-R@E]) if rows[other]['owner']=='elbow' or rows[other].get('nominal_pose_owner')=='elbow' else body;volume=intersection(candidate,target_body)
                if volume>.001:hits.append({'pair':[key,other],'volume_mm3':volume})
        records.append({'elbow_deg':angle,'relative_cap_swivel_deg':float(np.degrees(psi)),'body_axis_alignment_error_deg':float(np.degrees(np.arccos(np.clip((Q@z)@target,-1,1)))),'rotor_pose_affine_mm':T.tolist(),'rigid_intersections':hits})
    # Local fork is an actual material candidate, not a physical mounting claim.
    # Its next required connection is a real bone-rooted bridge; old integrated
    # cap/yoke cannot remain around this moving rotor.
    side=unit(np.cross(trunnion_axis,z));fork_trials=[];selected=None
    fork_parameters=itertools.product([-1,1],[16,20,26,32],[0,-8,8,-16,16,-24,24,-32,32]) if '--fork-offsets' in sys.argv else [(s,r,0) for s,r in itertools.product([-1,1],[20,26,32,40,48])]
    fork_parameters=sorted(fork_parameters,key=lambda t:(t[1]**2+t[2]**2,t))
    for sign,reach,depth in fork_parameters:
        direction=side*sign;ear_centres=[C+trunnion_axis*s*16.7 for s in [-1,1]];pieces=[]
        bridge_delta=direction*reach+z*depth
        for s,centre in zip([-1,1],ear_centres):
            pieces.append(rod(C+trunnion_axis*s*15.2,C+trunnion_axis*s*18.2,4))
            frame=np.column_stack([trunnion_axis,direction,np.cross(trunnion_axis,direction)])
            slab=mf.Manifold.cube([3,100,100]).translate([15.2 if s==1 else -18.2,-50,-50]).transform(np.c_[frame,C])
            # Keep the connecting web inside the same actual 3mm axial ear
            # thickness. A rounded web must not bulge into the bearing gap or
            # over the washer seating face.
            pieces.append(path([centre,centre+bridge_delta],2.5)^slab)
        pieces.append(rod(ear_centres[0]+bridge_delta,ear_centres[1]+bridge_delta,2.5));fork=mf.Manifold.batch_boolean(pieces,mf.OpType.Add)-rod(C-trunnion_axis*25,C+trunnion_axis*25,1.1)
        parts={'fork':fork}
        for s in [-1,1]:
            v=trunnion_axis*s
            pin_frame=np.column_stack([z,np.cross(v,z),v])
            # One stepped profile avoids coplanar cylinder-union slivers in
            # CAD sewing. Same nominal shaft/head lengths and radii.
            parts[f'pin_{s}']=bolt(7,18.7,1,2,1.3,32).transform(np.c_[pin_frame,C])
            parts[f'washer_{s}']=rod(C+v*18.2,C+v*18.7,2.5)-rod(C+v*17,C+v*20,1.1)
        static_pairs=[]
        for (ka,sa),(kb,sb) in itertools.combinations(parts.items(),2):
            v=intersection(sa,sb)
            if v>.001:static_pairs.append({'pair':[ka,kb],'volume_mm3':v})
        first=None
        for pose in records:
            R=Rotation.from_rotvec(axis*np.radians(pose['elbow_deg'])).as_matrix();Tf=np.c_[R,E-R@E];Tcap=np.array(pose['rotor_pose_affine_mm']);hits=[]
            for key,s in parts.items():
                current=s.transform(Tf)
                for other,body in moving.items():
                    v=intersection(current,body.transform(Tcap))
                    if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
                for other,body in fixed.items():
                    other_s=body.transform(Tf) if rows[other]['owner']=='elbow' or rows[other].get('nominal_pose_owner')=='elbow' else body;v=intersection(current,other_s)
                    if v>.001:hits.append({'pair':[key,other],'volume_mm3':v})
            if hits:first={'elbow_deg':pose['elbow_deg'],'hits':hits};break
        trial={'bridge_direction_display':direction.tolist(),'bridge_reach_mm':reach,'bridge_axial_shift_mm':depth,'bridge_delta_display_mm':bridge_delta.tolist(),'local_fixed_part_intersections':static_pairs,'first_conditional_interference':first};fork_trials.append(trial)
        if not static_pairs and first is None and all(len([c for c in s.decompose() if c.volume()>.001])==1 for s in parts.values()):
            selected=trial;cadparts=[meta]
            for key,s in parts.items():cadparts.append(write(dest/(key+'.stl'),s))
            (dest/'cad_input.json').write_text(json.dumps({'units':'mm','parts':cadparts},indent=2),encoding='utf8');break
    assert all(hashlib.sha256(__import__('pathlib').Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    report={'scope':'new physical rotor/fork and distal hardware vs existing rigid arm, conditional prescribed swivel; original integrated yoke absent and still requires an actual replacement bone-mounted bridge; NOT whole assembly clearance','source_sha256':hashes,'pivot_zero_display_mm':C.tolist(),'trunnion_axis_zero_display':trunnion_axis.tolist(),'stub_outer_radius_mm':1.8,'blind_pin_bore_radius_mm':1.1,'stub_axial_interval_mm':[10,15],'blind_bore_interval_mm':[11,16],'local_rotor_hardware_intersections':local_hits,'poses':records,'fork_trials':fork_trials,'selected_fork':selected,'missing':['bone-mounted support and its swept clearance','actual shoulder-pin threading/retention and bearing tolerances','TPU end deformation, wire/guide bending and neutral port fit','force balance, friction, preload, strength and physical testing'],'installed_mechanism_verified':False,'frozen':False};(dest/'screen.json').write_text(json.dumps(report,indent=2),encoding='utf8');print('swivel rotor localhits',len(local_hits),'conditional clear',sum(not r['rigid_intersections'] for r in records),'/',len(records),'forkselected',bool(selected),'forktrials',fork_trials,flush=True)

if __name__=='__main__':run()
