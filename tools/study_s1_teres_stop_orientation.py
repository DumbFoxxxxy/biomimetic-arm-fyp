"""Screen stop clocking using archived exact failed source geometry.

Coarse orientation search only. Selected orientation needs independent replay.
"""
from pathlib import Path
import numpy as np
from build_s1_teres_origin_mount import ROOT,A,rod,read,sha,save,load_component


def main():
    folder=A.DESIGN/'s1_teres_wire_stop';file=folder/'iterations/counterhold_blocked/assembly_check.json'
    r=read(file);I=np.array(r['local_origin_source_mm']);F=np.array(r['local_basis_columns']);n=F[:,2]
    solids={}
    for key,row in r['sources'].items():
        # The historical manifest retains its original paths. For its three
        # new parts load the byte-identical archived copies, so regenerating
        # the selected orientation does not invalidate this comparison.
        frozen=dict(row)
        if key in r['new_keys']:
            frozen['file']=(file.parent/Path(row['file']).name).relative_to(ROOT).as_posix()
        solids[key]=load_component(frozen)
    body_key='stop:FJ1507/hex_body';screw_key='stop:FJ1507/M2x3'
    fixed={key:value for key,value in solids.items() if key not in [body_key,screw_key]}
    mf=A.I._manifold();box=lambda lo,hi:mf.Manifold.cube(np.array(hi)-lo).translate(lo)
    fork=(box([-14.,-4.15,5.3],[-6.,4.15,6.8])+box([-8.,2.65,5.3],[1.,4.15,6.8])+
          box([-8.,-4.15,5.3],[1.,-2.65,6.8])).transform(np.c_[F,I])
    driver=rod(I+F@np.array([3.4,0.,6.05]),I+F@np.array([25.,0.,6.05]),.8)
    rows=[]
    angles=sorted(set(range(0,360,15))|set(range(175,196)))
    for angle in angles:
        rotation=A.I.rotation_about(I,n,float(angle))
        tools={'body':solids[body_key].transform(rotation),'screw':solids[screw_key].transform(rotation),
               'fork':fork.transform(rotation),'driver':driver.transform(rotation)}
        hits=[]
        for key,solid in tools.items():
            for target,obstacle in fixed.items():
                volume=float(A.I.intersection_volume(solid,obstacle))
                if volume>.001:hits.append(dict(part=key,target=target,volume_mm3=volume))
        path_hits=[]
        if not hits:
            direction=rotation[:,:3]@F[:,0]
            for distance in np.arange(0.,20.001,1.):
                tool=tools['fork'].translate(-direction*distance)
                for target,obstacle in fixed.items():
                    volume=float(A.I.intersection_volume(tool,obstacle))
                    if volume>.001:path_hits.append(dict(offset_mm=float(distance),target=target,volume_mm3=volume))
        rows.append(dict(angle_deg=angle,static_hits=hits,coarse_fork_path_hits=path_hits,
                         candidate_clear=not hits and not path_hits))
        print(angle,'static',len(hits),'path',len(path_hits),flush=True)
    save(folder/'orientation_study.json',dict(source_report=file.relative_to(ROOT).as_posix(),source_sha256=sha(file),
        script_sha256=sha(Path(__file__)),archived_new_meshes_resolved_beside_source_report=True,
        angle_step_deg=15,refined_angle_interval_deg=[175,195],refined_angle_step_deg=1,
        coarse_translation_step_mm=1.,rows=rows,
        not_verified='Independent final-source placement, real tools, grip, torque or drive'))


if __name__=='__main__':main()
