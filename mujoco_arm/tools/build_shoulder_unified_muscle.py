"""Monolithic TPU muscle/terminal candidate from the shared master wire paths.

Original anatomy, lattice and bone seats stay as reference. Dense terminal
sleeves and continuous channel walls are regenerated from the master paths.
This is static geometry, never a force, strain or contact-equilibrium model.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import mesh_geom
import interference as I
import check_shoulder_assembly as C
import build_shoulder_tendon_segments as T
import build_shoulder_ends as BE
import compact_shoulder_terminal as CT


def cap_keepout(inner, normal, thickness, radius):
    """Actual cap housing with an explicit 0.15mm fit gap on both faces.

    No arbitrary extension above the cap: it would remove the free TPU wall
    without there being a physical component to replace that wall.
    """
    inner=np.asarray(inner,float);normal=D.unit(normal)
    return T.rod(inner-normal*.15,inner+normal*(thickness+.15),radius+.15)


def discard_microscopic_positive_islands(solid):
    """Remove only minute detached positive debris; retain every cavity shell.

    The absolute and relative bounds are both required. Real disconnected
    pieces remain a failure; negative shells must never be filled as a repair.
    """
    pieces=solid.decompose()
    positive=sorted([p for p in pieces if p.volume()>0],key=lambda p:p.volume(),reverse=True)
    tiny=positive[1:]
    volume=sum(p.volume() for p in tiny)
    if not tiny or volume>=1e-3 or volume/max(1,solid.volume())>=1e-5:return solid,[]
    keep=[p for p in pieces if p.volume()<=0]+positive[:1]
    # compose()/union treats a negative-volume component as another solid
    # in this library version and fills cavities. Reassemble oriented mesh
    # shells without any Boolean operation instead.
    vertices=[];faces=[];offset=0
    for p in keep:
        v,f=mesh_geom.manifold_mesh_to_arrays(p.to_mesh64())
        vertices.append(v);faces.append(f+offset);offset+=len(v)
    rebuilt=I._manifold().Manifold(I._manifold().Mesh64(
        np.vstack(vertices),np.vstack(faces).astype(np.uint32)))
    removed=[dict(volume_mm3=p.volume(),bounds_local_mm=list(p.bounding_box())) for p in tiny]
    return rebuilt,removed


def local_export(solid,path,origin):
    local=I.clean_solid(solid.translate(-origin))
    original_volume=local.volume();removed=[]
    for attempt in range(3):
        v,f=mesh_geom.manifold_mesh_to_arrays(local.to_mesh())
        mesh_geom.write_stl(str(path),v,f,header='Unified TPU candidate, not tested performance')
        V,F=mesh_geom.read_stl(str(path));W,G=mesh_geom.weld(V,F)
        _,bnd,nm=mesh_geom.edge_census(G)
        restored=I._manifold().Manifold(I._manifold().Mesh64(W,G.astype(np.uint32)))
        cleaned,discarded=discard_microscopic_positive_islands(restored)
        if discarded and attempt<2:
            # A minute separate shell can itself cause a coincident edge in
            # the geometric census. Remove only the bounded recorded debris,
            # then check the written mesh again; never accept the bad census.
            local=cleaned;removed.extend(discarded)
            continue
        if (bnd or nm) and attempt<2:
            # Same bounded 0.001mm cleanup after the actual float32 export;
            # acceptance stays zero boundary and nonmanifold edges.
            local=I.clean_solid(restored,eps=1e-3)
            continue
        break
    components=[s.volume() for s in restored.decompose()]
    positive=[v for v in components if v>1e-9]
    cavities=[v for v in components if v< -1e-9]
    error=abs(restored.volume()-original_volume)/max(1,original_volume)
    return restored.translate(origin),dict(file=str(path.relative_to(D.frames.PROJECT_ROOT)),
        translation_source_mm=origin.tolist(),boundary_edges=bnd,nonmanifold_edges=nm,
        signed_component_volumes_mm3=components,roundtrip_relative_error=error,
        removed_positive_debris=removed,export_passes=attempt+1,
        positive_material_components=len(positive),enclosed_void_shells=len(cavities),
        mesh_checks_pass=(bnd==nm==0 and len(positive)==1 and error<1e-5),
        cavity_printability_verified=False,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--module',default='deltoid',
        choices=['deltoid','infraspinatus','subscapularis','teres_major'])
    ap.add_argument('--routes',type=Path,default=D.DEST/'routes_separated.json')
    ap.add_argument('--fit-terminals',action='store_true',help='Flat TPU mouths and remove classified nonstructural lattice cut-offs')
    args=ap.parse_args()
    source=args.routes;spec=json.loads(source.read_text(encoding='utf8'))
    routes=[r for r in spec['routes'] if r['module']==args.module]
    allseats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    seats=[s for s in allseats if s['module']==args.module]
    suffix='_outlet_socket' if any('rigid_origin_outlet' in r for r in routes) else '_fit_gap'
    if any('origin_stack_reduction_mm' in r for r in routes):suffix='_quarter_socket'
    if any('origin_seat_override' in r for r in routes):suffix='_upper_origin'
    edge_offsets=[r['origin_connector_edge_offset_mm'] for r in routes if 'origin_connector_edge_offset_mm' in r]
    if edge_offsets:suffix+='_edge'+str(edge_offsets[0]).replace('.','p')
    compact_routes=any('distal_terminal_override' in r for r in routes)
    if compact_routes:suffix='_compact_distal'
    if args.fit_terminals:suffix='_fitted'
    if compact_routes and args.fit_terminals:suffix='_compact_distal_fitted'
    folder=D.DEST/('unified_'+args.module+suffix);folder.mkdir(exist_ok=True)
    parent_file=D.frames.OUTPUT_DIR/'shoulder_tpu'/(args.module+'_tpu.stl')
    parent=C.load_stl(parent_file);pieces=[parent];caps={};keepouts=[];seams=[];assigned={};terminal_zones={};protected=[];stack_revisions={};attachment_revisions={};stops={};compact_specs={};bolt_drills=[]
    for s in seats:
        match=next(r for r in routes if r['head']==s['head'])
        override=match.get('origin_seat_override') if s['end']=='origin' else None
        compact=match.get('distal_terminal_override') if s['end']=='distal' else None
        if compact:s={**s,'cap_thickness_mm':compact['thickness_mm']}
        if override:
            old_seat=dict(s);s={**s,**override}
            attachment_revisions[s['head']+'/'+s['end']]=dict(original_centre_source_mm=old_seat['centre_source_mm'],
                original_normal_source=old_seat['normal_source'],replacement=override)
        filename=','.join(s['heads']).replace('/','_').replace(',','+')
        key=s['head']+'/'+s['end'];n=D.unit(s['normal_source'])
        assigned[key]=[h.split('/')[0] for h in s['heads']]
        port=match['ports'][s['end']]
        inner=np.array(port['inner_face_source_mm'])-np.array(port['offset_from_seat_mm'])
        # Prototype triplex cap/pad has true separate terminations. Other
        # source end parts remain referenced; their fastening is not validated.
        if len(s['heads'])==3:
            meta=json.loads((D.DEST/'triplex_terminal/check.json').read_text(encoding='utf8'))
            for kind in ['TPU_pad','rigid_cap']:
                row=meta['parts'][kind]
                body=C.load_stl(D.DEST/'triplex_terminal'/row['file']).translate(row['translation_source_mm'])
                if kind=='TPU_pad':pieces.append(body);protected.append(body)
                else:caps[key]=body
        else:
            if compact:
                cap,stop,specification=CT.build(inner,n,compact['tangent_source'])
                caps[key]=cap;stops[s['head']+'/distal_stop_blank']=stop;compact_specs[key]=specification
                c=np.array(s['centre_source_mm']);height=float((inner-c)@n)-.2
                if height<=0:raise ValueError('No TPU pad thickness remains')
                # Match the soft pad to the real two-lobe pressure plate;
                # retain both the steel and the closed M3 mounting regions.
                pad=CT.housing(c-n,n,compact['tangent_source'],radius=4.2,thickness=height+1)
                bolt=np.array(specification['bolt_inner_source_mm'])
                bolt_drills.append(T.rod(bolt-n*(height+3),bolt+n*6,1.7))
            elif override:
                c=np.array(s['centre_source_mm']);basis=D.shoulder_pads.pad_basis(n)
                pad=D.shoulder_pads.build_pad(c-n,n,thickness=s['pad_thickness_mm']+1,radius=s['pad_radius_mm'])
                for sign in [-1,1]:
                    p=c+basis[:,0]*7*sign
                    pad=pad-T.rod(p-n*2,p+n*(s['pad_thickness_mm']+3),1.7)
                caps[key]=BE.disc(inner,n,s['cap_thickness_mm'],radius=s['cap_radius_mm'],holes=True)
            else:
                pad=C.load_stl(D.frames.OUTPUT_DIR/'shoulder_parts'/(filename+'_pad.stl'))
                caps[key]=C.load_stl(D.frames.OUTPUT_DIR/'shoulder_parts'/(filename+'_cap.stl'))
            reduction_key=s['end']+'_stack_reduction_mm'
            if reduction_key in match:
                reduction=float(match[reduction_key]);maximum=3.+(s['lift_mm'] if s['end']=='distal' else 0.)
                if not 0<=reduction<=maximum+1e-9:raise ValueError('Stack reduction exceeds boss and lift')
                if not compact:caps[key]=caps[key].translate(-n*reduction)
                plane=inner-n*.2
                fitted_pad=pad.trim_by_plane(-n,float((-n)@plane))
                stack_revisions[key]=dict(reduction_mm=reduction,
                    pad_material_removed_mm3=pad.volume()-fitted_pad.volume(),
                    pad_outer_plane_source_mm=plane.tolist(),bone_seat_moved=False)
                pad=fitted_pad
            pieces.append(pad);protected.append(pad)
        # Trim soft overgrowth locally above the clamping face. This is a
        # defined cap footprint, not an entire-bone collision exclusion.
        mask=CT.housing(inner-n*.15,n,compact['tangent_source'],radius=3.35,thickness=3.3) if compact else cap_keepout(inner,n,s['cap_thickness_mm'],s['cap_radius_mm'])
        keepouts.append(mask)
        seams.append(T.rod(inner-n*.15,inner,s['cap_radius_mm']+.15))
        terminal_zones[key]=mask
    outer_guards=[];lumens=[];wire={};guard={};lumens_by_head={};outlet_caps=[];outlet_specs={}
    for r in routes:
        h=r['head'];full=np.array(r['full_wire_path_source_mm'])
        lumen=D.outer_envelope(full,.86);lumens.append(lumen)
        lumens_by_head[h]=lumen
        wire[h]=D.outer_envelope(full,.3)
        if outlet:=r.get('rigid_origin_outlet'):
            key=h+'/origin';outlet_path=np.array(outlet['path_source_mm'])
            actual_radius=float(outlet['outer_radius_mm'])
            caps[key]=caps[key]+D.outer_envelope(outlet_path,actual_radius)
            # Physical fit candidate, not a tolerance on intersection results.
            socket=D.outer_envelope(outlet_path,actual_radius+.10)
            keepouts.append(socket);terminal_zones[key]=terminal_zones[key]+socket
            outlet_caps.append(key);outlet_specs[key]=dict(**outlet,socket_radial_gap_mm=.10)
        core=np.array(r['core_path_source_mm'])
        # Core OD3.6mm gives nominal wall0.94mm after the 1.72mm lumen.
        core_wall=D.outer_envelope(core,1.8);pieces.append(core_wall);protected.append(core_wall)
        core_guard=D.outer_envelope(core,1.61)-lumen
        endpoint_guards=[]
        for end in ['origin','distal']:
            path=np.array(r['connectors'][end]['path_source_mm'])
            sleeve=D.outer_envelope(path,3.2)
            if args.fit_terminals:
                # Rounded envelope endpoints can leave a free TPU nub on the
                # far side of a thin cap. Terminate this sleeve at the actual
                # inner interface plane instead of trying to reconnect that nub.
                n=D.unit(r['ports'][end]['normal_source'])
                sleeve=sleeve-T.rod(path[0]-n*.001,path[0]+n*6.4,3.21)
            pieces.append(sleeve);protected.append(sleeve)
            # The first4mm is the separately checked cap/terminal zone, not
            # the free compliant tube. Check its actual wire clearance below.
            length=np.linalg.norm(np.diff(path,axis=0),axis=1).sum()
            if length>4.1:
                region=D.trim_path(path,4.,length)
                endpoint_guards.append(D.outer_envelope(region,1.61)-lumen)
        guard[h]=I._manifold().Manifold.batch_boolean([core_guard]+endpoint_guards,I._manifold().OpType.Add)
    print(args.module,'union lattice/pads/walls/transitions',flush=True)
    solid=I._manifold().Manifold.batch_boolean(pieces,I._manifold().OpType.Add)
    bones={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()}
    bone_union=I._manifold().Manifold.batch_boolean(list(bones.values()),I._manifold().OpType.Add)
    trim_volume=float((solid^bone_union).volume());solid=solid-bone_union
    cap_bone_trims={key:I.intersection_volume(caps[key],bone_union) for key in outlet_caps}
    for key in outlet_caps:caps[key]=caps[key]-bone_union
    for mask in keepouts:solid=solid-mask
    hole_union=I._manifold().Manifold.batch_boolean(lumens,I._manifold().OpType.Add)
    solid=solid-hole_union
    for drill in bolt_drills:solid=solid-drill
    classified_tips=[]
    if args.fit_terminals:
        functional=I._manifold().Manifold.batch_boolean(protected,I._manifold().OpType.Add)
        components=sorted([p for p in solid.decompose() if p.volume()>1e-9],key=lambda p:p.volume(),reverse=True)
        for fragment in components[1:]:
            volume=fragment.volume()
            protected_hit=I.intersection_volume(fragment,functional)
            lattice_fraction=I.intersection_volume(fragment,parent)/volume
            # No size-only deletion gate. Only detached pieces from the old
            # lattice that touch NONE of the pads, active cores or sleeves.
            if protected_hit<=1e-5 and lattice_fraction>=.9999:
                classified_tips.append(dict(volume_mm3=volume,bounds_source_mm=list(fragment.bounding_box()),
                    protected_material_intersection_mm3=protected_hit,source_lattice_fraction=lattice_fraction,
                    reason='detached original lattice cut-off; no attachment, active channel or transition'))
                solid=solid-fragment
    # Every own cap is drilled from the same routes. Preserve source parts.
    caps={k:s-I._manifold().Manifold.batch_boolean(
        [lumens_by_head[h] for h in assigned[k]],I._manifold().OpType.Add) for k,s in caps.items()}
    origin=np.array(seats[0]['centre_source_mm'])
    solid,soft_meta=local_export(solid,folder/'TPU_muscle.stl',origin)
    parts={'TPU':solid};metadata={'TPU':soft_meta}
    for key,s in caps.items():
        parts[key],metadata[key]=local_export(s,folder/(key.replace('/','_')+'_cap.stl'),origin)
    for key,s in stops.items():
        parts[key],metadata[key]=local_export(s,folder/(key.replace('/','_')+'.stl'),origin)
        metadata[key]['kind']='stop_blank';metadata[key]['attachment_body']='humerus'
    checks={h:{key:I.intersection_volume(w,s) for key,s in parts.items()} for h,w in wire.items()}
    wire_bones={h:{key:I.intersection_volume(w,s) for key,s in bones.items()} for h,w in wire.items()}
    # An internal clear lumen does not demonstrate an accessible external
    # mouth. Check the nominal hand-pull endpoint against actual materials.
    lead_endpoints={}
    for row in routes:
        tip=np.array(row['ports']['origin']['hand_pull_lead_source_mm'])
        probe=I._manifold().Manifold.sphere(.3,32).translate(tip)
        lead_endpoints[row['head']]={key:I.intersection_volume(probe,s)
            for key,s in {**parts,**{'bone:'+k:v for k,v in bones.items()}}.items()}
    material_pairs=[dict(a=a,b=b,intersection_mm3=I.intersection_volume(parts[a],parts[b]))
        for a,b in itertools.combinations(parts,2)]
    wall={h:float((g-solid).volume()) for h,g in guard.items()}
    wall_free={}
    for h,g in guard.items():
        zones=[z for k,z in terminal_zones.items() if h in assigned[k]]
        excluded=I._manifold().Manifold.batch_boolean(zones,I._manifold().OpType.Add)
        missing=(g-solid)-excluded
        wall_free[h]=dict(missing_mm3=float(missing.volume()),
                         bounds_source_mm=list(missing.bounding_box()) if missing.volume()>1e-6 else None)
    # Other muscles are not silently cut to make the new design fit.
    neighbours={}
    for module in ['deltoid','infraspinatus','subscapularis','teres_major']:
        if module!=args.module:
            body=C.load_stl(D.frames.OUTPUT_DIR/'shoulder_tpu'/(module+'_tpu.stl'))
            neighbours[module]=I.intersection_volume(solid,body)
    report=dict(scope='static monolithic TPU candidate from master paths; not physical drive or stage acceptance',
        module=args.module,master_routes=str(source.relative_to(D.frames.PROJECT_ROOT)),
        master_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(parent_file.read_bytes()).hexdigest(),
        parts=metadata,wire_material_intersections_mm3=checks,missing_0_75mm_wall_guard_mm3=wall,
        wire_bone_intersections_mm3=wire_bones,
        nominal_hand_pull_endpoint_intersections_mm3=lead_endpoints,
        material_pair_intersections_mm3=material_pairs,
        cap_keepout_scope='actual cap thickness plus explicit 0.15mm fit gap at both faces and radius; actual rigid outlet socket where present',
        rigid_origin_outlets=outlet_specs,rigid_cap_bone_trim_mm3=cap_bone_trims,
        origin_stack_revisions=stack_revisions,
        attachment_revisions=attachment_revisions,
        compact_distal_terminals=compact_specs,
        compact_stop_capture_inward_mm3={key:{str(d):I.intersection_volume(parts[key+'/distal_stop_blank'].translate(
            -D.unit(spec['normal_source'])*d),parts[key+'/distal']) for d in [0,.04,.06]}
            for key,spec in [(k.split('/')[0],v) for k,v in compact_specs.items()]},
        compact_bone_mounting_holes_implemented=False,
        fitted_terminal_mouths=args.fit_terminals,removed_nonstructural_lattice_tips=classified_tips,
        all_part_bone_intersections_mm3={key:{k:I.intersection_volume(part,b) for k,b in bones.items()} for key,part in parts.items()},
        free_TPU_wall_guard_outside_assigned_cap_housings=wall_free,
        required_wall_mm=.75,wall_guard_scope='retained active core and connector arclength>=4mm; cap fit housing and actual rigid guide socket excluded explicitly',
        dense_core_outer_radius_mm=1.8,transition_outer_radius_mm=3.2,lumen_radius_mm=.86,
        conforming_bone_trim_volume_mm3=trim_volume,bone_intersections_mm3={k:I.intersection_volume(solid,s) for k,s in bones.items()},
        unchanged_neighbour_intersections_mm3=neighbours,stage_pass=False,
        not_verified=['external lead accessibility and bone clearance','neighbour resolution and complete assembly','fasteners/crimps and installation','material stroke, friction and deformation','joint retention and motion'])
    (folder/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Mesh checks',{k:r['mesh_checks_pass'] for k,r in metadata.items()},flush=True)
    print('Wire hits',checks,'wall guard missing',wall,'free TPU guard',wall_free,'neighbours',neighbours,flush=True)
    print('Wire to bone',wire_bones,'lead endpoints',lead_endpoints,flush=True)


if __name__=='__main__':main()
