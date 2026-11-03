"""Explicit local opening experiment for isolated non-channel voids.

Never fill negative shells or mutate S1-Y1. Nominal aperture sizes are
unvalidated geometry trials, not a printer specification or strength pass.
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

DEST=A.DESIGN/'s1_y_muscle_pair/vent_experiment'

def positive_void(shell):
    v,f=mesh_geom.manifold_mesh_to_arrays(shell.to_mesh64())
    return A.I._manifold().Manifold(A.I._manifold().Mesh64(
        np.array(v,dtype=np.float64,order='C',copy=True),
        np.array(f[:,::-1],dtype=np.uint64,order='C',copy=True)))

def negative_volume(solid):
    return -sum(p.volume() for p in solid.decompose() if p.volume() < -1e-6)

def search(current,void,protected,radius):
    v,f=mesh_geom.manifold_mesh_to_arrays(void.to_mesh64());tri=v[f]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm=np.linalg.norm(cross,axis=1)
    normals=cross/np.maximum(norm[:,None],1e-30);centres=tri.mean(axis=1)
    # Face-area ordered, deterministic spatially distinct seeds. Search at
    # most 80 local directions; failure means this search found no opening.
    seeds=[]
    for i in np.argsort(-norm):
        if not any(np.linalg.norm(centres[i]-centres[j])<radius*.5 for j in seeds):seeds.append(int(i))
        if len(seeds)>=80:break
    rows=[];best=None;old_negative=negative_volume(current)
    for i in seeds:
        q=centres[i];n=normals[i]
        hits=current.ray_cast(q+n*.02,q+n*12.)
        ds=sorted(float((np.array(h.position)-q)@n) for h in hits if float((np.array(h.position)-q)@n)>.04)
        if not ds:continue
        length=ds[0];start=q-n*radius;end=q+n*(length+radius)
        cutter=rod(start,end,radius);removed=current^cutter
        item=dict(face=i,radius_mm=radius,path_source_mm=[start.tolist(),end.tolist()],
            wall_distance_mm=length,removed_material_mm3=float(removed.volume()),
            protected_material_removed_mm3=A.I.intersection_volume(removed,protected))
        if item['protected_material_removed_mm3']>1e-6:
            item['rejection']='enters protected endpoint or R1.8 internal channel stock';rows.append(item);continue
        trial=current-cutter;pieces=trial.decompose()
        item['positive_components']=sum(p.volume()>1e-6 for p in pieces)
        item['remaining_enclosed_volume_mm3']=negative_volume(trial)
        if item['positive_components']!=1:item['rejection']='disconnects material'
        elif old_negative-item['remaining_enclosed_volume_mm3']<void.volume()*.9:item['rejection']='does not open intended cavity'
        else:
            item['accepted_for_geometry_comparison']=True
            if best is None or item['removed_material_mm3']<best[1]['removed_material_mm3']:best=(trial,item)
        rows.append(item)
    return best,rows

def run():
    file=A.DESIGN/'s1_y_muscle_pair/assembly_check.json';assembly=json.loads(file.read_text(encoding='utf8'))
    routefile=A.DESIGN/'s1_port_path_revision/routes_compact_distal_candidate.json'
    routes=json.loads(routefile.read_text(encoding='utf8'))['routes'];results=[]
    DEST.mkdir(parents=True,exist_ok=True)
    for module in ['infraspinatus','subscapularis']:
        source=assembly['sources']['tpu:'+module];old,_,_=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])
        route=next(r for r in routes if r['module']==module)
        protection=[]
        for p in route['ports'].values():
            c=np.array(p['inner_face_source_mm']);n=np.array(p['normal_source']);protection.append(rod(c-n*8.,c+n*3.,12.))
        for path in [route['core_path_source_mm']]+[route['connectors'][side]['path_source_mm'] for side in ['origin','distal']]:
            protection.append(A.D.outer_envelope(np.array(path),1.8))
        protected=old^A.I._manifold().Manifold.batch_boolean(protection,A.I._manifold().OpType.Add)
        original_voids=[positive_void(s) for s in old.decompose() if s.volume() < -1e-6]
        # Each aperture size starts independently from the unchanged source.
        for radius in [.4,.6,.8]:
            current=old;steps=[]
            for i,void in enumerate(sorted(original_voids,key=lambda v:-v.volume())):
                best,attempts=search(current,void,protected,radius)
                step=dict(cavity_index=i,cavity_volume_mm3=float(void.volume()),attempts=attempts)
                if best is None:step['result']='no acceptable local opening in finite search'
                else:current,step['selected']=best;step['result']='local opening found'
                steps.append(step)
            folder=DEST/module/('diameter_'+f'{2*radius:g}');folder.mkdir(parents=True,exist_ok=True)
            origin=np.array(source['translation_source_mm'])
            local,mesh=write_readback(current.translate(-origin),folder/'TPU_vent_trial.stl')
            actual=local.translate(origin)
            mesh['translation_source_mm']=origin.tolist()
            mesh['mesh_checks_pass']=(mesh['boundary_edges']==mesh['nonmanifold_edges']==0
                and mesh['positive_components']==1 and mesh['relative_volume_roundtrip_error']<1e-5)
            result=dict(module=module,aperture_diameter_candidate_mm=2*radius,parent_TPU=source,
                steps=steps,mesh=mesh,parent_closed_void_volume_mm3=negative_volume(old),
                remaining_closed_void_volume_mm3=negative_volume(actual),
                net_removed_material_mm3=float((old-actual).volume()),numerical_added_material_mm3=float((actual-old).volume()),
                protected_material_lost_mm3=float((protected-actual).volume()),
                all_cavities_opened_in_search=all('selected' in s for s in steps),
                manufacturing_strength_or_drive_verified=False,stage_pass=False)
            results.append(result)
            print(module,'diameter',2*radius,'remaining void',result['remaining_closed_void_volume_mm3'],
                  'mesh',mesh['mesh_checks_pass'],'opened',result['all_cavities_opened_in_search'],flush=True)
    report=dict(scope=__doc__,assembly_source=file.relative_to(A.ROOT).as_posix(),assembly_sha256=A.hashed(file),
        routes_source=routefile.relative_to(A.ROOT).as_posix(),routes_sha256=A.hashed(routefile),
        script_sha256=A.hashed(Path(__file__).resolve()),
        export_helper_sha256=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_Y_print_inspection.py'),
        results=results,adopted_as_current=False,
        stage_pass=False,actual_drive_range_deg=None)
    (DEST/'vent_trials.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
