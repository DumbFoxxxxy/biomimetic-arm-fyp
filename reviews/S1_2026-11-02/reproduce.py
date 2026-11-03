"""Read-only S1 audit. Run with Python environment providing numpy/manifold3d."""
from pathlib import Path
import json
import sys
import hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'.deps'),str(ROOT/'mujoco_arm/tools')]
import frames
import interference as I
import build_shoulder_ends as BE
import check_shoulder_assembly as C
import sweep_gh_bare as B


def main():
    p=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'
    d=json.loads((p/'s1_modules.json').read_text(encoding='utf8'))
    recorded=json.loads((p/'s1_cable_length.json').read_text(encoding='utf8'))
    G=frames.joint_centres('source')['GH'];axes=frames.gh_probe_axes()
    result={'scope':'independent audit, no production geometry modified','force_lines':[],'pads':[],'bone_at_120':[]}
    for axis,head in [('flexion','FJ1467'),('abduction','FJ1468'),('longitudinal','FJ1500')]:
        c=next(c for m in d['modules'] for c in m['corridors'] if c['head']==head)
        A=np.array(c['path_source_mm'][0]);b=np.array(c['path_source_mm'][-1]);k=axes[axis]
        ts=np.arange(0,120.0001,.1)
        moved=np.array([(lambda T:T[:,:3]@b+T[:,3])(I.rotation_about(G,k,t)) for t in ts])
        L=np.linalg.norm(moved-A,axis=1)
        # Positive means tension produces positive rotation on the probed axis.
        arm=-np.einsum('ij,ij->i',(moved-A)/L[:,None],np.cross(k,moved-G))
        cross=np.where(arm[:-1]*arm[1:]<0)[0]
        old=next(r for r in recorded['rows'] if r['axis']==axis and r['head']==head)
        result['force_lines'].append(dict(axis=axis,head=head,recorded_monotonic=old['monotonic'],
            sign_change_brackets_deg=[[float(ts[i]),float(ts[i+1])] for i in cross],
            moment_arm_zero_mm=float(arm[0]),moment_arm_120_mm=float(arm[-1]),
            moment_arm_min_mm=float(arm.min()),moment_arm_max_mm=float(arm.max()),
            chord_start_mm=float(L[0]),chord_end_mm=float(L[-1])))
    folder=frames.OUTPUT_DIR/'shoulder_parts'
    for f in sorted(folder.glob('*pad.stl')):
        s=C.load_stl(f);pieces=[a for a in s.decompose() if a.volume()>.001]
        result['pads'].append(dict(file=str(f.relative_to(ROOT)),components=len(pieces),
            volumes_mm3=[float(a.volume()) for a in pieces],sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    primitive=BE.disc([0,0,0],[0,0,1],4,boss_mm=3)
    result['pad_primitive_bounds']=[list(a.bounding_box()) for a in primitive.decompose()]
    moving=B.bone_solid(B.HUMERUS)
    fixed=[B.bone_solid(B.SCAPULA),B.bone_solid(B.CLAVICLE)]
    for axis in ['abduction','longitudinal']:
        s=moving.transform(I.rotation_about(G,axes[axis],120))
        result['bone_at_120'].append(dict(axis=axis,scapula_overlap_mm3=I.intersection_volume(s,fixed[0]),
            clavicle_overlap_mm3=I.intersection_volume(s,fixed[1])))
    # Preserve results.json as the initial audit snapshot after production fixes.
    (Path(__file__).parent/'results_current.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print('Audit results saved; disconnected pads:',sum(r['components']!=1 for r in result['pads']))


if __name__=='__main__':main()
