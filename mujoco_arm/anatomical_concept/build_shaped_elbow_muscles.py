"""Anatomy-derived PSM manufacturing candidates: lattice/shell/guides/end pads.
Engineering contraction coordinates are not verified anatomical fiber directions.
"""
from parameterize_elbow_muscles import *
from build_SC_connection import rod,write
from scipy.spatial import cKDTree

PARAM={'strut_diameter_mm':2.,'cell_spacing_transverse_mm':[15.,20.],'cell_spacing_axial_mm':10.,'surface_strands_each_direction':8,'channel_radius_mm':.8,'guide_outer_radius_mm':2.,'guide_length_mm':3.,'guide_spacing_mm':10.,'end_pad_radius_mm':11.,'end_pad_thickness_mm':4.,'end_bolt_spacing_mm':14.,'end_bolt_hole_radius_mm':1.7,'wire_diameter_candidate_mm':.6}

def capsule(a,b,r):return rod(a,b,r)+mf.Manifold.sphere(r,12).translate(a)+mf.Manifold.sphere(r,12).translate(b)

def run(muscles=('FJ1486','FJ1477')):
    for fid in muscles:
        profile=json.loads((DEST/(fid+'_profile.json')).read_text(encoding='utf8'));folder=DEST/fid;folder.mkdir(exist_ok=True)
        secs=profile['sections'];u=np.array([s['height_mm'] for s in secs]);contours=np.array([s['boundary_local_mm'] for s in secs]);poles=np.array([s['pole_local_mm'] for s in secs]);F=np.array(profile['basis_columns_display']);E=np.array(profile['coordinate_origin_display_mm'])
        curve=make_interp_spline(u,contours,k=3,axis=0);poles=route_engineering_poles(secs)
        def centre(h):return np.array([np.interp(h,u,poles[:,i]) for i in [0,1]])
        def ring(h,n=64):
            xy=curve(h);xy=np.vstack([xy,xy[0]])
            return make_interp_spline(np.linspace(0,1,33),xy,k=3,bc_type='periodic')(np.linspace(0,1,n,endpoint=False))
        def surface(h,t):
            xy=curve(h);xy=np.vstack([xy,xy[0]]);p=make_interp_spline(np.linspace(0,1,33),xy,k=3,bc_type='periodic')(t%1)
            return np.r_[p,h]
        def c(h):return np.r_[centre(np.clip(h,u[0],u[-1])),h]
        heights=np.linspace(u[0],u[-1],int(np.ceil((u[-1]-u[0])/3))+1);n=64
        vertices=np.array([np.c_[ring(h),np.full(n,h)] for h in heights]).reshape(-1,3);faces=[]
        for i in range(len(heights)-1):
            for j in range(n):
                a=i*n+j;b=i*n+(j+1)%n;d=(i+1)*n+j;e=(i+1)*n+(j+1)%n
                faces.extend([[a,b,e],[a,e,d]])
        faces.extend(np.asarray(mf.triangulate([ring(heights[0])]))[:,::-1].tolist())
        faces.extend((np.asarray(mf.triangulate([ring(heights[-1])]))+(len(heights)-1)*n).tolist())
        body=mf.Manifold(mf.Mesh64(vertices.astype(np.float64),np.array(faces,dtype=np.uint64)))
        if body.status()!=mf.Error.NoError:raise RuntimeError(fid+' invalid sectional loft '+str(body.status()))
        if body.volume()<0:raise RuntimeError('Negative loft orientation')
        # Cartesian eight-branch cuboid cells, trimmed by the actual curved loft.
        spacing=np.array([15.,20.,10.]);lo=np.floor(vertices.min(0)/spacing).astype(int);hi=np.ceil(vertices.max(0)/spacing).astype(int);rods=[]
        for indices in itertools.product(*(range(a,b) for a,b in zip(lo,hi))):
            corner=np.array(indices)*spacing;mid=corner+spacing/2
            for bits in itertools.product([0,1],repeat=3):rods.append(rod(mid,corner+np.array(bits)*spacing,1.))
        internal=mf.Manifold.batch_boolean(rods,mf.OpType.Add)^body
        # A controlled surface diamond net, unrelated to source mesh triangles.
        surface_rods=[];layers=np.linspace(u[0],u[-1],int(np.ceil((u[-1]-u[0])/5))+1)
        for sign in [-1,1]:
            for strand in range(8):
                points=[surface(h,strand/8+sign*i/16) for i,h in enumerate(layers)]
                for a,b in zip(points,points[1:]):
                    steps=max(2,int(np.ceil(np.linalg.norm(b-a)/2))+1)
                    # Sample the parametric curve rather than a long contour chord.
                    index=int(np.searchsorted(layers,a[2]));h0,h1=a[2],b[2];t0=strand/8+sign*index/16
                    samples=[surface(h0+(h1-h0)*f,t0+sign*f/16) for f in np.linspace(0,1,steps)]
                    surface_rods.extend(rod(p,q,1.) for p,q in zip(samples,samples[1:]))
        shell=mf.Manifold.batch_boolean(surface_rods,mf.OpType.Add)
        stations=np.linspace(u[0],u[-1],int(np.ceil((u[-1]-u[0])/10))+1);guide_parts=[];spokes=[];guide_specs=[]
        for index,h in enumerate(stations):
            p=c(h);previous=c(stations[index-1]) if index else p-[0,0,9];following=c(stations[index+1]) if index+1<len(stations) else p+[0,0,9]
            incoming=unit(p-previous);outgoing=unit(following-p);tangent=unit(incoming+outgoing);alpha=np.arccos(np.clip(tangent@incoming,-1,1))
            radial_cut=.85/max(.01,np.cos(alpha));outer=max(2.,radial_cut+.6+.4*np.tan(alpha)+.1)
            length=min(3.,2*(outer-radial_cut-.6)/max(1e-10,np.tan(alpha)))
            if length<.8:raise RuntimeError('Excessive nominal wire turn '+fid)
            guide_parts.append(rod(p-tangent*length/2,p+tangent*length/2,outer))
            guide_specs.append({'centre_display_mm':(p@F.T+E).tolist(),'axis_display':(F@tangent).tolist(),'length_mm':float(length),'outer_radius_mm':float(outer),'incoming_outgoing_turn_deg':float(np.degrees(2*alpha)),'nominal_tube_wall_lower_bound_mm':float(outer-radial_cut-length/2*np.tan(alpha)),'interpretation':'conservative nominal straight-span envelope only; deformation/friction not verified'})
            for t in [0,.25,.5,.75]:spokes.append(rod(p,surface(h,t),1.))
        guides=mf.Manifold.batch_boolean(guide_parts,mf.OpType.Add)
        webs=mf.Manifold.batch_boolean(spokes,mf.OpType.Add)^body
        muscle=internal+shell+guides+webs
        # Dense TPU pads and branched transitions are actual material, with M3 holes.
        pads=[];pad_centres=[];transitions=[]
        for sign,h in [(-1,u[0]),(1,u[-1])]:
            p=c(h)+[0,0,sign*9];pad_centres.append(p)
            pads.append(rod(p-[0,0,2],p+[0,0,2],11.))
            for t in np.linspace(0,1,16,endpoint=False):transitions.append(rod(p+[7*np.cos(2*np.pi*t),7*np.sin(2*np.pi*t),-sign*1.5],surface(h,t),1.2))
        muscle=muscle+mf.Manifold.batch_boolean(pads+transitions,mf.OpType.Add)
        hardware={};nut_seats=[];nut_access_windows=[]
        for label,sign,p in [('proximal',-1,pad_centres[0]),('distal',1,pad_centres[1])]:
            for side,x in [('left',-7),('right',7)]:
                q=p+[x,0,0];head=q+[0,0,sign*6.65]
                hardware[label+'_bolt_'+side]=rod(head,head-[0,0,sign*12],1.5)+rod(head,head+[0,0,sign*3],2.75)
                hardware[label+'_head_washer_'+side]=rod(q+[0,0,sign*6.15],head,3.5)-rod(q-[0,0,16],q+[0,0,16],1.7)
                hardware[label+'_nut_washer_'+side]=rod(q-[0,0,sign*2],q-[0,0,sign*2.5],3.5)-rod(q-[0,0,16],q+[0,0,16],1.7)
                frame=np.diag([1.,-sign,-sign]);base=q-[0,0,sign*2.5]
                nut=mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).transform(np.c_[frame,base])-rod(q-[0,0,16],q+[0,0,16],1.6)
                hardware[label+'_nut_'+side]=nut
                nut_seats.append(mf.Manifold.batch_boolean([nut]+[nut.translate(d*s*.15) for d in np.eye(3) for s in [-1,1]],mf.OpType.Add))
                washer=hardware[label+'_nut_washer_'+side]
                nut_seats.append(mf.Manifold.batch_boolean([washer]+[washer.translate(d*s*.15) for d in np.eye(3) for s in [-1,1]],mf.OpType.Add))
                # Actual outward-facing access slot, not merely a final-state
                # recess: admit the nut+washer before the bolt is installed.
                direction=np.sign(x);zlo,zhi=sorted([-sign*5.1,-sign*1.8])
                lo=np.array([min(x,x+direction*20)-3.8,-3.8,zlo]);hi=np.array([max(x,x+direction*20)+3.8,3.8,zhi])
                nut_seats.append(mf.Manifold.cube(hi-lo).translate(p+lo))
                nut_access_windows.append({'end':label,'side':side,'bounds_local_mm':np.r_[p+lo,p+hi].tolist(),'insertion_direction_local':[-float(direction),0,0],'distance_mm':20,'pad_material_policy':'slot begins 0.2mm below inner pad face; remaining nominal pad thickness3.8mm, no strength claim'})
        terminal_clearance=mf.Manifold.batch_boolean(nut_seats,mf.OpType.Add)
        muscle=muscle-terminal_clearance
        # Close any generated disconnected branches by real local material bridges.
        fragments=sorted(muscle.decompose(),key=lambda s:s.volume(),reverse=True);bridges=[];main=fragments[0]
        if len(fragments)>100:raise RuntimeError('Too many disconnected branches '+fid)
        for fragment in fragments[1:]:
            if abs(fragment.volume())<=1e-10:continue
            av=np.asarray(main.to_mesh64().vert_properties)[:,:3];bv=np.asarray(fragment.to_mesh64().vert_properties)[:,:3];distance,index=cKDTree(av).query(bv);i=int(np.argmin(distance));a=bv[i];b=av[index[i]]
            if distance[i]>10:raise RuntimeError('Nonlocal lattice fragment '+fid)
            main=main+fragment+capsule(a,b,1.);bridges.append({'from_local_mm':a.tolist(),'to_local_mm':b.tolist(),'gap_mm':float(distance[i])})
        muscle=main
        # Wire is straight between actual short guide tubes. No invisible spline constraint.
        wire_nodes=[pad_centres[0]-[0,0,20],*[c(h) for h in stations],pad_centres[1]+[0,0,8]]
        nominal_channel=mf.Manifold.batch_boolean([capsule(a,b,.8) for a,b in zip(wire_nodes,wire_nodes[1:])],mf.OpType.Add)
        channel=mf.Manifold.batch_boolean([capsule(a,b,.85) for a,b in zip(wire_nodes,wire_nodes[1:])],mf.OpType.Add)
        holes=[]
        for p in pad_centres:
            for sign in [-1,1]:holes.append(rod(p+[sign*7,0,-8],p+[sign*7,0,8],1.7))
        bolt_void=mf.Manifold.batch_boolean(holes,mf.OpType.Add)
        muscle=muscle-channel-bolt_void-terminal_clearance
        # Drilling can detach an otherwise connected generated branch. Reconnect
        # with explicit reinforced geometry while retaining BOTH designed voids.
        post_cut_bridges=[];fragments=sorted(muscle.decompose(),key=lambda s:s.volume(),reverse=True);main=fragments[0]
        for fragment in fragments[1:]:
            if abs(fragment.volume())<=1e-10:continue
            av=np.asarray(main.to_mesh64().vert_properties)[:,:3];bv=np.asarray(fragment.to_mesh64().vert_properties)[:,:3];distance,index=cKDTree(av).query(bv);connected=False
            for i in np.argsort(distance)[:24]:
                if distance[i]>10:break
                a=bv[i];b=av[index[i]]
                bridge=capsule(a,b,2.)-channel-bolt_void-terminal_clearance;trial=main+fragment+bridge
                positive=[s for s in trial.decompose() if abs(s.volume())>1e-10]
                if len(positive)==1:
                    main=positive[0];connected=True;post_cut_bridges.append({'from_local_mm':a.tolist(),'to_local_mm':b.tolist(),'radius_mm':2.,'gap_mm':float(distance[i])});break
            if not connected:raise RuntimeError('Unable to reconnect drilled branch '+fid)
        muscle=main
        muscle=muscle.as_original().simplify(.01)
        pore_fills=[];sliver_trims=[]
        for component in muscle.decompose():
            volume=float(component.volume());bounds=np.array(component.bounding_box());span=np.linalg.norm(bounds[3:]-bounds[:3])
            if -0.01<volume< -1e-10:
                mesh=component.to_mesh64();fill=mf.Manifold(mf.Mesh64(np.asarray(mesh.vert_properties)[:,:3].astype(np.float64),np.asarray(mesh.tri_verts)[:,::-1].astype(np.uint64)))
                muscle=muscle+fill;pore_fills.append({'filled_generated_closed_pore_mm3':-volume})
            elif 1e-10<volume<.001 and span<1.5:
                centre_trim=(bounds[:3]+bounds[3:])/2;radius_trim=span/2+.05
                muscle=muscle-mf.Manifold.sphere(radius_trim,16).translate(centre_trim);sliver_trims.append({'centre_local_mm':centre_trim.tolist(),'radius_mm':radius_trim,'generated_sliver_mm3':volume})
        components=[s for s in muscle.decompose() if s.volume()>1e-6]
        # Cap files contain real cable ports; bone adapters and metal fasteners remain separate work.
        caps={}
        for label,sign,p in [('proximal',-1,pad_centres[0]),('distal',1,pad_centres[1])]:
            a=p+[0,0,sign*2.15];b=a+[0,0,sign*4];cap=rod(a,b,11.)
            for x in [-7,7]:cap=cap-rod(p+[x,0,-12],p+[x,0,12],1.7)
            cap=cap-rod(p-[0,0,16],p+[0,0,16],.8)
            if label=='proximal':cap=cap+rod(b,b-[0,0,6],3.)-rod(b+[0,0,1],b-[0,0,7],1.1)
            else:
                cap=cap-rod(b-[0,0,3.2],b+[0,0,.1],1.4)
                cap=cap-rod((a+b)/2+[-11,4,0],(a+b)/2+[11,4,0],1.)
                hardware['active_wire_crimp_barrel']=rod(b-[0,0,3.05],b-[0,0,.05],1.25)-rod(b-[0,0,4],b+[0,0,1],.45)
            caps[label+'_rigid_cap']=cap.as_original().simplify(.005)
        parts={'TPU_muscle':muscle,**caps,**hardware};rows=[]
        for name,s in parts.items():
            if s.status()!=mf.Error.NoError:raise RuntimeError('Invalid '+fid+' '+name)
            transformed=s.transform(np.c_[F,E]);rows.append(write(folder/(name+'.stl'),transformed))
        sourcev,sourcef=read_obj(str(ROOT/profile['source']));source=mf.Manifold(mf.Mesh64(((sourcev-SHIFT-E)@F).astype(np.float64),sourcef.astype(np.uint64)))
        report={'id':fid,'scope':'shaped TPU manufacturing geometry candidate, NOT completed installed muscle/drive validation','parameters':PARAM,'profile_sha256':hashlib.sha256((DEST/(fid+'_profile.json')).read_bytes()).hexdigest(),'coordinate_origin_display_mm':E.tolist(),'basis_columns_display':F.tolist(),'body_span_mm':float(u[-1]-u[0]),'loft_volume_mm3':float(body.volume()),'loft_outside_original_anatomy_mm3':max(0,float((body-source).volume())),'TPU_material_volume_mm3':float(muscle.volume()),'material_inside_loft_fraction':float((muscle^body).volume()/body.volume()),'TPU_positive_component_count':len(components),'guide_centres_display_mm':[(c(h)@F.T+E).tolist() for h in stations],'wire_polyline_display_mm':[(p@F.T+E).tolist() for p in wire_nodes],'pad_centres_display_mm':[(p@F.T+E).tolist() for p in pad_centres],'real_branch_connectors':bridges,'components':{r['name']:r['components'] for r in rows},'active_force_path':'motor/Bowden reaction at proximal rigid cap -> active steel inside short guide tubes -> crimp stop in distal rigid cap; distal passive tendon through transverse eyelet connects to joint. TPU pads require actual M3 hardware and bone adapters.','pending':['M3 fastener stack and bone adapters','actual crimp/ferrule and conduit seats verification','joint tendon routing/moment arms and required stroke','nonlinear deformation, bone contact, guide bending and friction','printing, cyclic recovery, damage and load tests'],'fiber_direction_verified':False,'contraction_capability_verified':False,'frozen':False}
        report['post_cut_reinforcement_bridges']=post_cut_bridges
        report['generated_closed_pore_fills']=pore_fills;report['generated_sliver_trims']=sliver_trims
        report['fasteners']='Four nominal M3x12 bolts, four M3 nuts, eight washers; thread/preload not modeled. Metal crimp barrel must actually be crimped to grip wire.'
        report['nut_access_windows']=nut_access_windows
        report['actual_short_tube_guides']=guide_specs
        report['engineering_guide_poles_local_mm']=poles.tolist()
        report['wire_path_method']='fixed-end minimum-step path through interior section candidates, linear interpolation; no anatomical fiber claim or invisible guide constraints'
        report['pending']=['proximal rigid cap to bone adapter and motor/conduit reaction fixture','passive tendon route and distal bone/rotor attachment','actual crimp grip, fastener preload and TPU end creep','joint moment arms, required stroke and compliant deformation','bone/soft contact, guide friction, printing and cyclic damage tests']
        report['neutral_module_overlaps']=[]
        for (an,a),(bn,b) in itertools.combinations(parts.items(),2):
            volume=max(0,float((a^b).volume()))
            if volume>.001:report['neutral_module_overlaps'].append({'parts':[an,bn],'volume_mm3':volume})
        report['channel_cut_radius_mm']=.85;report['nominal_channel_minimum_target_radius_mm']=.8
        report['wire_channel_TPU_overlap_mm3']=max(0,float((muscle^nominal_channel).volume()))
        (folder/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');(folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':rows}),encoding='utf8')
        print(fid,'TPU components',len(components),'bridges',len(bridges),'volume',round(muscle.volume(),1),'caps',[r['components'] for r in rows[1:]],flush=True)

if __name__=='__main__':run()
