"""Independent Y-centreline lattice bellies on the current anatomical S1 model.

Keeps existing end-interface material and continuous channel/transition stock.
Rebuilds the remaining lattice in recorded transported paper XYZ frames.
No source paper performance, print result or physical muscle stroke is assumed.
"""
import json
import shutil
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import build_shoulder_tpu as T
from build_shoulder_unified_muscle import local_export
from build_shoulder_shaft_collar import rod
from check_shoulder_front_liner_keepers import is_compliant
from scipy.spatial import cKDTree

DEST=A.DESIGN/'s1_y_muscle_pair'
BASELINE=A.DESIGN/'s1_GH_assembly_review/compact_recessed_keepers/assembly_check.json'
RECIPES={'infraspinatus':[(15.,10.,30.),(10.,8.,10.)],
         'subscapularis':[(15.,10.,30.),(6.,8.,20.),(10.,8.,10.)]}


def connect_replacement_stock(solid,old,required,junctions,allowed,lumen,rigid):
    """Classify clipped tails by role, then bridge retained functional stock.

    No volume-only deletion rule: all required pad/hook/channel and new cell
    junction material remains. Added bridges have two positive node overlaps.
    """
    removed=[];bridges=[]
    for part in sorted(solid.decompose(),key=lambda p:-p.volume())[1:]:
        if part.volume()<=1e-9:continue
        volume=float(part.volume());functional=A.I.intersection_volume(part,required)
        old_fraction=A.I.intersection_volume(part,old)/volume
        node_hit=A.I.intersection_volume(part,junctions)
        if functional<=1e-5 and (old_fraction>=.999 or node_hit<=1e-5):
            removed.append(dict(volume_mm3=volume,old_material_fraction=old_fraction,
                                required_material_intersection_mm3=functional,new_junction_intersection_mm3=node_hit,
                                bounds_source_mm=list(part.bounding_box()),
                                reason='obsolete old-grid clipping remnant' if old_fraction>=.999 else 'new clipped strut tail containing no cell centre/corner or functional stock'))
            solid-=part
    for _ in range(30):
        pieces=sorted([p for p in solid.decompose() if p.volume()>1e-6],key=lambda p:-p.volume())
        if len(pieces)<=1:return solid,removed,bridges,None
        main=pieces[0];fragment=pieces[1]
        V=np.array(main.to_mesh64().vert_properties)[:,:3];W=np.array(fragment.to_mesh64().vert_properties)[:,:3]
        # Prefer preserving a bounded original interface transition to inventing
        # a new arch around a clamp. The stock has already fitted this hardware.
        centre=W.mean(0);chosen=None
        for radius in [3.,5.,8.,11.]:
            patch=(old^A.I._manifold().Manifold.sphere(radius,32).translate(centre))-lumen
            contacts=[A.I.intersection_volume(patch,s) for s in [main,fragment]]
            if min(contacts)<.02:continue
            trial=solid+patch
            if sum(x.volume()>1e-6 for x in trial.decompose())>=len(pieces):continue
            chosen=(trial,dict(method='bounded restoration of existing interface-transition stock',
                               source_centre_mm=centre.tolist(),source_radius_mm=radius,
                               fragment_volume_mm3=float(fragment.volume()),node_overlap_mm3=contacts,
                               added_material_mm3=float((trial-solid).volume())));break
        if chosen is not None:
            solid,row=chosen;bridges.append(row);continue
        W=np.vstack([W.mean(0),W[np.linspace(0,len(W)-1,min(100,len(W)),dtype=int)]])
        distance,index=cKDTree(V).query(W)
        chosen=None
        for j in np.argsort(distance)[:80]:
            p,q=W[j],V[index[j]];length=float(np.linalg.norm(p-q))
            if length>8.:continue
            bridge=T.tube([p,q],1.,segments=24)-lumen
            if float((bridge-allowed).volume())>.001:continue
            if any(A.I.intersection_volume(bridge,s)>.001 for s in rigid.values()):continue
            contacts=[A.I.intersection_volume(bridge,s) for s in [main,fragment]]
            if min(contacts)<.02:continue
            trial=solid+bridge
            if sum(x.volume()>1e-6 for x in trial.decompose())>=len(pieces):continue
            chosen=(trial,dict(method='new full-radius short transition strut',path_source_mm=[p.tolist(),q.tolist()],radius_mm=1.,length_mm=length,
                               fragment_volume_mm3=float(fragment.volume()),node_overlap_mm3=contacts,
                               added_material_mm3=float((trial-solid).volume())));break
        if chosen is None:return solid,removed,bridges,'no full-radius short bridge preserves anatomy/rigid/wire constraints'
        solid,row=chosen;bridges.append(row)
    return solid,removed,bridges,'component reconnection iteration limit'


def export_with_local_junction_blend(candidate,path,origin,lumen,rigid):
    """Repair only measured edge-only junctions using explicit local material.

    Never accept a nonmanifold mesh or compose negative shells to fill bores.
    Every attempted blend remains a real geometric edit, checked against wire
    lumen and rigid parts, and its location/radius/material delta is recorded.
    """
    actual,meta=local_export(candidate,path,origin);attempts=[]
    if meta['mesh_checks_pass']:return actual,meta,attempts
    if meta['boundary_edges'] or meta['nonmanifold_edges']>8 or meta['positive_material_components']!=1:
        return actual,meta,attempts
    V,F=A.mesh_geom.read_stl(str(path));V,F=A.mesh_geom.weld(V,F)
    edges=np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);edges.sort(axis=1)
    unique,count=np.unique(edges,axis=0,return_counts=True)
    centres=[V[p].mean(0)+origin for p in unique[count>2]]
    if not centres:return actual,meta,attempts
    archive=path.parent/'iterations/pinched_readback';archive.mkdir(parents=True,exist_ok=True)
    shutil.copy2(path,archive/path.name);(archive/'mesh.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf8')
    for radius in [.15,.3,.6,1.]:
        beads=T.union_all([A.I._manifold().Manifold.sphere(radius,32).translate(c) for c in centres])
        blend=beads-lumen
        if any(A.I.intersection_volume(blend,s)>.001 for s in rigid.values()):
            attempts.append(dict(radius_mm=radius,rejected='rigid material conflict'));continue
        trial=(candidate+blend)-lumen
        restored,row=local_export(trial,path,origin)
        attempts.append(dict(radius_mm=radius,centres_source_mm=[c.tolist() for c in centres],
                             added_material_mm3=float((trial-candidate).volume()),mesh_checks_pass=row['mesh_checks_pass'],
                             boundary_edges=row['boundary_edges'],nonmanifold_edges=row['nonmanifold_edges']))
        actual,meta=restored,row
        if row['mesh_checks_pass']:break
    return actual,meta,attempts


def thickness_seed(tri,path,envelope):
    P=np.asarray(path);length=np.linalg.norm(np.diff(P,axis=0),axis=1);s=np.r_[0,np.cumsum(length)]
    half=s[-1]/2.;i=min(max(0,np.searchsorted(s,half)-1),len(P)-2)
    t=(P[i+1]-P[i])/length[i];c=P[i]+t*(half-s[i]);V=tri.reshape(-1,3)
    selected=V[abs((V-c)@t)<5.]
    if len(selected)<10:raise ValueError('Insufficient middle-belly section vertices')
    x=A.D.unit(np.cross(np.eye(3)[np.argmin(abs(t))],t));z=np.cross(x,t)
    uv=np.column_stack([(selected-c)@x,(selected-c)@z]);values,basis=np.linalg.eigh(np.cov(uv.T))
    thin=A.D.unit(x*basis[0,0]+z*basis[1,0]);wide=np.cross(thin,t)
    def spans(axis):return sorted(float((np.array(h.position)-c)@axis) for h in envelope.ray_cast(c-axis*80,c+axis*80))
    return thin,dict(method='PCA of envelope mesh vertices within +/-5mm of middle arclength plane; geometric estimate only',
                     centre_source_mm=c.tolist(),section_normal_source=t.tolist(),vertex_count=len(selected),
                     covariance_eigenvalues_mm2=values.tolist(),paper_X_seed_source=thin.tolist(),
                     estimated_middle_X_intersections_mm=spans(thin),estimated_middle_Z_intersections_mm=spans(wide))


def run():
    DEST.mkdir(parents=True,exist_ok=True)
    selection_path=A.ROOT/'design/S1_CURRENT_CANDIDATE.json';selection=json.loads(selection_path.read_text(encoding='utf8'))
    # Fixed Z-belly reference keeps regeneration idempotent if this candidate
    # later becomes the selected assembly; do not rebuild on our own output.
    source=BASELINE
    parent=json.loads(source.read_text(encoding='utf8'))
    S={k:A.load_checked(A.ROOT/r['file'],r['translation_source_mm'],r['sha256'])[0] for k,r in parent['sources'].items()}
    route_file=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(route_file.read_text(encoding='utf8'))['routes'];sources=dict(parent['sources']);attempts=[];chosen={}
    for module,recipes in RECIPES.items():
        key='tpu:'+module;old=S[key];route=next(r for r in routes if r['module']==module)
        env_file=A.frames.OUTPUT_DIR/'shoulder_modules'/f'{module}.stl';envelope,tri,_=A.load_checked(env_file)
        path=np.array(route['core_path_source_mm']);seed,pca=thickness_seed(tri,path,envelope)
        ends=[];inner_ends=[]
        for side,port in route['ports'].items():
            c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source'])
            ends.append(rod(c-n*8.,c+n*3.,12.))
            inner_ends.append(rod(c-n*7.5,c+n*2.5,11.5))
        end_guard=T.union_all(ends);new_lattice_end_exclusion=T.union_all(inner_ends)
        full=np.array(route['full_wire_path_source_mm'])
        internal_paths=[path]+[np.array(route['connectors'][side]['path_source_mm']) for side in ['origin','distal']]
        channel_guard=T.union_all([A.D.outer_envelope(p,3.3) for p in internal_paths])
        lumen=A.D.outer_envelope(full,.86)
        protected=end_guard+channel_guard
        retained=old-(envelope-protected)
        origin_port=route['ports']['origin'];c=np.array(origin_port['inner_face_source_mm']);n=np.array(origin_port['normal_source'])
        pad_height=2.5 if route['head']=='FJ1500' else 4.
        functional=rod(c-n*pad_height,c,11.)
        port=route['ports']['distal'];c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source'])
        functional+=rod(c-n*4.,c,4.5)
        tangent=np.array(route['distal_terminal_override']['tangent_source']);tangent=A.D.unit(tangent-n*(tangent@n))
        F=np.column_stack([tangent,np.cross(n,tangent),n])
        for lo,hi in [([-1.,3.4,-2.],[1.,5.,4.]),([-1.,2.6,3.2],[1.,5.,4.]),
                      ([-1.,-5.,-2.],[1.,-3.4,4.]),([-1.,-5.,3.2],[1.,-2.6,4.])]:
            functional+=A.I._manifold().Manifold.cube((np.array(hi)-lo).tolist()).translate(lo).transform(np.c_[F,c])
        required_internal_channel=T.union_all([A.D.outer_envelope(p,1.8) for p in internal_paths])
        functional=(functional+required_internal_channel)^old
        for cell in recipes:
            folder=DEST/module/('_'.join(f'{v:g}' for v in cell));folder.mkdir(parents=True,exist_ok=True)
            frame={};lattice,cells,skipped=T.structured_lattice(envelope,tri,path,cell=cell,r=1.,path_axis='y',
                                                             initial_transverse_axis=seed,frame_record=frame)
            row=dict(module=module,head=route['head'],cell_xyz_mm=list(cell),strut_diameter_mm=2.,
                     working_axis='paper_y_along_recorded_centreline',pca_section=pca,transported_frames=frame,
                     envelope_file=env_file.relative_to(A.ROOT).as_posix(),envelope_sha256=A.hashed(env_file),
                     current_TPU_sha256=parent['sources'][key]['sha256'],cells=cells,skipped_cells=skipped,
                     parent_volume_mm3=float(old.volume()),retained_old_material_mm3=float(retained.volume()),
                     protected_end_cylinders=dict(radius_mm=12.,normal_span_from_inner_face_mm=[-8.,3.]),
                     transition_overlap_band_mm=.5,
                     new_lattice_end_exclusion=dict(radius_mm=11.5,normal_span_from_inner_face_mm=[-7.5,2.5]),
                     protected_internal_channel_radius_mm=3.3,
                     channel_guard_scope='retained core and origin/distal connector paths; external hand-pull lead excluded',
                     selected=False,stage_pass=False,actual_drive_range_deg=None)
            if lattice is None:
                row['rejection']='no structured cells fit this anatomical envelope at these dimensions'
                attempts.append(row);(folder/'design.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf8')
                print(module,cell,'REJECT no cells',flush=True);continue
            # A finite overlap band joins two distinct lattice grids. A zero
            # overlap cut can leave an edge-only contact and a pinched cavity.
            new_stock=(lattice^envelope)-new_lattice_end_exclusion
            removals=[]
            for obstacle_key,obstacle in S.items():
                if is_compliant(obstacle_key):continue
                volume=A.I.intersection_volume(new_stock,obstacle)
                if volume>1e-9:
                    before=new_stock
                    # A local nominal relief is a geometry choice, not print tolerance.
                    lo,hi=np.array(new_stock.bounding_box()).reshape(2,3)
                    crop=A.I._manifold().Manifold.cube((hi-lo+.4).tolist()).translate(lo-.2)
                    keepout=(obstacle^crop).minkowski_sum(A.I._manifold().Manifold.cube([.2]*3,True))
                    new_stock-=keepout
                    removals.append(dict(part=obstacle_key,initial_overlap_mm3=volume,removed_new_stock_mm3=float((before-new_stock).volume())))
            raw=retained+new_stock;candidate=raw-lumen
            nodes=frame['cell_centres_source_mm']+frame['live_corner_nodes_source_mm']
            junctions=T.union_all([A.I._manifold().Manifold.sphere(.35,12).translate(p) for p in nodes])
            candidate,obsolete,bridges,reconnection_failure=connect_replacement_stock(
                candidate,old,functional,junctions,envelope+end_guard+channel_guard,lumen,
                {k:s for k,s in S.items() if not is_compliant(k)})
            positive=sorted([float(p.volume()) for p in candidate.decompose() if p.volume()>1e-6],reverse=True)
            row.update(rigid_reliefs=removals,new_lattice_volume_mm3=float(new_stock.volume()),
                       raw_positive_component_volumes_mm3=positive,removed_by_lumen_mm3=float((raw-(raw-lumen)).volume()),
                       removed_classified_tails=obsolete,interface_transition_bridges=bridges,
                       reconnection_failure=reconnection_failure)
            actual,meta,blend_attempts=export_with_local_junction_blend(candidate,folder/'TPU_Y_muscle.stl',
                np.array(parent['sources'][key]['translation_source_mm']),lumen,
                {k:s for k,s in S.items() if not is_compliant(k)})
            row['mesh']=meta
            row['local_junction_blends']=blend_attempts
            row['wire_overlap_mm3']=A.I.intersection_volume(A.D.outer_envelope(full,.3),actual)
            row['protected_end_material_lost_mm3']=float(((old^end_guard)-actual).volume())
            row['protected_channel_material_lost_mm3']=float(((old^channel_guard)-actual).volume())
            row['required_pad_hook_channel_material_lost_mm3']=float((functional-actual).volume())
            row['retained_old_fraction_of_parent']=float(retained.volume()/old.volume())
            row['actual_new_material_mm3']=float((actual-old).volume());row['actual_removed_material_mm3']=float((old-actual).volume())
            row['rigid_hits']=[dict(part=k,volume_mm3=A.I.intersection_volume(actual,s)) for k,s in S.items() if not is_compliant(k)]
            row['rigid_hits']=[r for r in row['rigid_hits'] if r['volume_mm3']>.001]
            row['accepted_geometry_candidate']=bool(meta['mesh_checks_pass'] and not row['rigid_hits'] and row['wire_overlap_mm3']<=.001)
            if row['accepted_geometry_candidate']:
                row['selected']=True;S[key]=actual;sources[key]=meta;chosen[module]=row
            else:row['rejection']='disconnected/readback-invalid lattice, rigid overlap, or blocked nominal wire; retain diagnostic geometry'
            attempts.append(row);(folder/'design.json').write_text(json.dumps(row,indent=2)+'\n',encoding='utf8')
            print(module,cell,'cells',cells,'components',len(positive),'mesh',meta['mesh_checks_pass'],
                  'rigid hits',row['rigid_hits'],'accepted',row['accepted_geometry_candidate'],flush=True)
            if row['selected']:break
    report=dict(scope=__doc__,parent_report=source.relative_to(A.ROOT).as_posix(),parent_sha256=A.hashed(source),
                selection_sha256=A.hashed(selection_path),script_sha256=A.hashed(Path(__file__).resolve()),routes_sha256=A.hashed(route_file),
                generation_dependency_sha256={p.relative_to(A.ROOT).as_posix():A.hashed(p) for p in
                    [A.ROOT/'mujoco_arm/tools/build_shoulder_tpu.py',A.ROOT/'mujoco_arm/tools/build_shoulder_unified_muscle.py',A.ROOT/'mujoco_arm/tools/build_shoulder_shaft_collar.py']},
                attempts=attempts,selected_modules=list(chosen),both_muscles_have_geometry_candidate=len(chosen)==2,
                sources=sources,n_instances=len(sources),stage_pass=False,actual_drive_range_deg=None,
                not_verified=['independent all-pair readback, interface/channel guard and frame checks',
                              'surface protective pattern completeness, print orientation and slicing',
                              'continuous-channel axial stiffness, deformable contact and free-terminal/tendon equilibrium',
                              'real shortening, force, friction, recovery and physical drive'])
    (DEST/'candidate_generation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Selected muscles',list(chosen),flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');run()
