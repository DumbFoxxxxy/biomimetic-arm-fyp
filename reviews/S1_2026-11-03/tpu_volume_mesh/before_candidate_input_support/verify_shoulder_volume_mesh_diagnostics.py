"""Independent algebra/topology checks of rejected TPU tetrahedralization input."""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A

BASE=A.DESIGN/'s1_y2_transition_stock/volume_mesh'

def run():
    file=BASE/'diagnostics.json';d=json.loads(file.read_text(encoding='utf8'))
    if d['script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/diagnose_shoulder_tpu_volume_mesh.py'):raise ValueError('Changed diagnostic code')
    if d['mesher_script_sha256']!=A.hashed(A.ROOT/'mujoco_arm/tools/build_shoulder_tpu_volume_mesh.py'):raise ValueError('Changed mesher')
    cube=np.load(BASE/'control_cube.npz');nodes=cube['nodes_local_mm'];tet=cube['tetrahedra'];p=nodes[tet]
    volumes=np.linalg.det(np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=2))/6.
    if np.any(volumes<=0) or abs(volumes.sum()-1000.)>1e-6:raise ValueError('Independent cube volume failed')
    faces=np.vstack([tet[:,[0,1,2]],tet[:,[0,1,3]],tet[:,[0,2,3]],tet[:,[1,2,3]]])
    unique,counts=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True)
    if counts.max()>2 or np.sum(counts==1)!=12:raise ValueError('Control cube tetrahedral conformity failed')
    intersections=[];sources=[]
    for row in d['localized_failures']:
        folder=(A.ROOT/row['open_patch']).parent
        if A.hashed(folder/'engine_stdout.txt')!=row['native_stdout_sha256']:raise ValueError('Changed native engine log')
        attempt=json.loads((folder/'mesh_attempt.json').read_text(encoding='utf8'));source=attempt['source']
        if A.hashed(A.ROOT/source['file'])!=source['sha256'] or not attempt['source_bytes_unchanged']:raise ValueError('Source changed')
        if attempt['geometry_tetrahedralization_pass'] or attempt['numerical_elasticity_mesh_quality_pass']:raise ValueError('Failure promoted to a pass')
        for crossing in row['strict_nonadjacent_segment_facet_crossings']:
            p,q=np.array(crossing['edge_points_local_mm']);a,b,c=np.array(crossing['facet_points_local_mm'])
            matrix=np.column_stack([q-p,-(b-a),-(c-a)]);values=np.linalg.lstsq(matrix,a-p,rcond=None)[0]
            t,u,v=values;residual=float(np.linalg.norm(matrix@values-(a-p)));condition=float(np.linalg.cond(matrix))
            point=p+t*(q-p);difference=float(np.linalg.norm(point-crossing['intersection']['point_local_mm']))
            if residual>1e-9 or difference>1e-8 or not (-1e-7<=t<=1+1e-7 and u>1e-8 and v>1e-8 and u+v<1-1e-8):
                raise ValueError('Independent local intersection check failed')
            intersections.append(dict(module=row['module'],least_squares_parameters=values.tolist(),residual_mm=residual,
                matrix_condition=condition,point_difference_mm=difference))
        sources.append(dict(module=row['module'],source_sha256=source['sha256'],
            raw_exact_nonmanifold_edges=attempt['input_preprocessing']['raw_nonmanifold_edges'],
            imported_nonmanifold_edges=attempt['input_preprocessing']['preprocessed_nonmanifold_edges'],
            source_coordinate_failure_point_mm=row['source_coordinate_point_mm'],
            engine_point_inside_Y2_infill_bounds=row['engine_point_inside_any_Y2_infill_bounds']))
    for attempt in d['attempts']:
        if A.hashed(A.ROOT/attempt['file'])!=attempt['sha256']:raise ValueError('Attempt report changed')
        if attempt['geometry_pass'] or attempt['elasticity_quality_pass']:raise ValueError('Unexpected complete-muscle pass')
    report=dict(scope=__doc__,diagnostics_sha256=A.hashed(file),script_sha256=A.hashed(Path(__file__).resolve()),
        control_cube=dict(tetrahedra=len(tet),positive_volume_mm3=float(volumes.sum()),boundary_faces=int(np.sum(counts==1)),passed=True),
        independently_checked_local_intersections=intersections,source_mesh_checks=sources,
        all_attempts_remain_rejected=True,complete_muscle_elasticity_ready=False,
        stage_pass=False,actual_drive_range_deg=None)
    (BASE/'verification_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print('Cube control',len(tet),'tets; independently checked intersections',len(intersections),'; full-muscle mesh gate remains failed.',flush=True)

if __name__=='__main__':run()
