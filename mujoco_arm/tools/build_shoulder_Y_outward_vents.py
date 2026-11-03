"""Finite vent search with cavity-bounded backward penetration.

The prior cutter extended one aperture radius behind each cavity face. That
can cross a thin pocket's opposite wall. Here the start depth is derived
from the actual positive empty-volume ray; protected stock is unchanged.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A
import mesh_geom
from build_shoulder_shaft_collar import rod
from build_shoulder_Y_print_inspection import write_readback
from build_shoulder_Y_vent_candidates import positive_void,negative_volume

DEST=A.DESIGN/'s1_y_muscle_pair/vent_experiment/outward_start_revision'

def search(current,void,protected,radius):
    v,f=mesh_geom.manifold_mesh_to_arrays(void.to_mesh64());tri=v[f]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);sizes=np.linalg.norm(cross,axis=1)
    normals=cross/np.maximum(sizes[:,None],1e-30);centres=tri.mean(axis=1)
    seeds=[]
    for i in np.argsort(-sizes):
        if not any(np.linalg.norm(centres[i]-centres[j])<radius*.5 for j in seeds):seeds.append(int(i))
        if len(seeds)>=80:break
    attempts=[];best=None;previous=negative_volume(current)
    for i in seeds:
        q=centres[i];n=normals[i]
        empty_hits=void.ray_cast(q-n*1e-5,q-n*12.)
        widths=sorted(float((q-np.array(h.position))@n) for h in empty_hits if float((q-np.array(h.position))@n)>1e-4)
        if not widths:continue
        thickness=widths[0];back=min(.02,thickness*.2)
        hits=current.ray_cast(q+n*1e-5,q+n*12.)
        lengths=sorted(float((np.array(h.position)-q)@n) for h in hits if float((np.array(h.position)-q)@n)>.0001)
        if not lengths:continue
        length=lengths[0];start=q-n*back;end=q+n*(length+radius)
        cutter=rod(start,end,radius);removed=current^cutter
        prior_cutter=rod(q-n*radius,end,radius)
        item=dict(face=i,radius_mm=radius,face_point_source_mm=q.tolist(),outward_normal_source=n.tolist(),
            cavity_thickness_along_normal_mm=thickness,backward_start_depth_mm=back,
            previous_backward_depth_mm=radius,path_source_mm=[start.tolist(),end.tolist()],wall_distance_mm=length,
            removed_material_mm3=float(removed.volume()),protected_material_removed_mm3=A.I.intersection_volume(removed,protected),
            previous_same_face_protected_loss_mm3=A.I.intersection_volume(current^prior_cutter,protected))
        if item['protected_material_removed_mm3']>1e-6:
            item['rejection']='protected stock intersection';attempts.append(item);continue
        trial=current-cutter;item['positive_components']=sum(p.volume()>1e-6 for p in trial.decompose())
        item['remaining_enclosed_volume_mm3']=negative_volume(trial)
        if item['positive_components']!=1:item['rejection']='disconnects material'
        elif previous-item['remaining_enclosed_volume_mm3']<void.volume()*.9:item['rejection']='does not open intended pocket'
        else:
            item['accepted_for_geometry_comparison']=True
            if best is None or item['removed_material_mm3']<best[1]['removed_material_mm3']:best=(trial,item)
        attempts.append(item)
    return best,attempts

def run():
    file=A.DESIGN/'s1_y_muscle_pair/assembly_check.json';assembly=json.loads(file.read_text(encoding='utf8'))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json';routes=json.loads(routefile.read_text(encoding='utf8'))['routes']
    source=assembly['sources']['tpu:subscapularis'];old,_,_=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])
    r=next(r for r in routes if r['module']=='subscapularis');regions=[]
    for port in r['ports'].values():
        c=np.array(port['inner_face_source_mm']);n=np.array(port['normal_source']);regions.append(rod(c-n*8.,c+n*3.,12.))
    for path in [r['core_path_source_mm']]+[r['connectors'][s]['path_source_mm'] for s in ['origin','distal']]:
        regions.append(A.D.outer_envelope(np.array(path),1.8))
    protected=old^A.I._manifold().Manifold.batch_boolean(regions,A.I._manifold().OpType.Add)
    voids=sorted([positive_void(s) for s in old.decompose() if s.volume() < -1e-6],key=lambda p:-p.volume())
    results=[];DEST.mkdir(parents=True,exist_ok=True)
    for radius in [.4,.6,.8]:
        actual=old;steps=[]
        for i,void in enumerate(voids):
            best,attempts=search(actual,void,protected,radius)
            step=dict(cavity_index=i,cavity_volume_mm3=float(void.volume()),attempts=attempts)
            if best is None:step['result']='no acceptable opening in finite search'
            else:actual,step['selected']=best;step['result']='opening found with cavity-bounded start'
            steps.append(step)
        folder=DEST/('diameter_'+f'{radius*2:g}');folder.mkdir(exist_ok=True)
        origin=np.array(source['translation_source_mm']);local,mesh=write_readback(actual.translate(-origin),folder/'TPU_vent_trial.stl')
        restored=local.translate(origin);mesh['translation_source_mm']=origin.tolist()
        result=dict(module='subscapularis',aperture_diameter_candidate_mm=2*radius,parent_TPU=source,
            steps=steps,mesh=mesh,parent_closed_void_volume_mm3=negative_volume(old),
            remaining_closed_void_volume_mm3=negative_volume(restored),
            net_removed_material_signed_mm3=float((old-restored).volume()),
            numerical_added_material_signed_mm3=float((restored-old).volume()),
            protected_material_lost_signed_mm3=float((protected-restored).volume()),
            all_cavities_opened_in_search=all('selected' in s for s in steps),stage_pass=False)
        results.append(result)
        print('subscapularis diameter',2*radius,'remaining void',result['remaining_closed_void_volume_mm3'],
              'opened',result['all_cavities_opened_in_search'],flush=True)
    report=dict(scope=__doc__,assembly_source=file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(file),
        routes_source=routefile.relative_to(A.ROOT).as_posix(),routes_sha256=A.hashed(routefile),
        script_sha256=A.hashed(Path(__file__).resolve()),
        dependency_sha256={name:A.hashed(A.ROOT/'mujoco_arm/tools'/name) for name in ['build_shoulder_Y_vent_candidates.py','build_shoulder_Y_print_inspection.py']},
        results=results,adopted_as_current=False,stage_pass=False,actual_drive_range_deg=None)
    (DEST/'outward_trials.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
