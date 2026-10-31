"""Index elbow design candidates honestly; do not freeze untested actuation."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'mujoco_arm/exports/anatomical_concept'
BASE=OUT/'direct_muscle_interfaces/mechanical_yoke_reanchored_pair/paired_swivel_support/captive_keepers'
CONT=BASE/'continuous_channel_candidate'


def run():
    a=json.loads((CONT/'assembly.json').read_text())
    rows={r['key']:r for r in a['instances']}
    interfaces=[]
    for fid in ['FJ1486','FJ1477']:
        d=json.loads((CONT/fid/'design.json').read_text());offset=np.array(rows[fid+'_TPU_muscle']['translation_display_mm'])
        A,B=np.array(d['pad_centres_display_mm'])+offset;F=np.array(d['basis_columns_display']);z=F[:,2]
        interfaces.append(dict(id=fid,frame='display millimetres; raw BP3D minus [-140,-70,1100]',
            proximal_TPU_pad_centre_display_mm=A.tolist(),distal_TPU_pad_centre_display_mm=B.tolist(),
            proximal_conduit_reaction_shoulder_display_mm=(A-z*5.15).tolist(),
            distal_crimp_nominal_seat_display_mm=(B+z*4.6).tolist(),
            wire_threading_reference_polyline_display_mm=(np.array(d['wire_polyline_display_mm'])+offset).tolist(),
            wire_threading_note='Original polyline includes threading overrun beyond crimp seat; trim working terminal at nominal crimp seat, not final polyline point.',
            distal_pivot_nominal_centre_display_mm=(B+z*4.15).tolist(),
            proximal_owner='humerus', distal_owner=fid+'_distal_pivot',
            actual_distal_pivot_instances=[r['key'] for r in a['instances'] if r['owner']==fid+'_distal_pivot'],
            basis_columns_display=F.tolist(),assembly_translation_display_mm=offset.tolist(),
            nominal_pose_only=True,passive_equilibrium_verified=False,crimp_grip_verified=False))
    (OUT/'elbow_interfaces.json').write_text(json.dumps(dict(scope='authoritative nominal interface metadata from actual selected sources; supersedes inherited prose muscle_end_boundaries',interfaces=interfaces),indent=2),encoding='utf8')
    inputs=[CONT/'assembly.json',CONT/'candidate_review.json',OUT/'elbow_test_coupons/manifest.json',OUT/'elbow_interfaces.json',BASE/'review_status.json']
    for fid in ['FJ1486','FJ1477']: inputs += [CONT/fid/'TPU_muscle.obj',CONT/fid/'TPU_muscle.stl']
    evidence=dict(source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        actual_TPU_drive_range_degrees=None,physical_performance_verified=False,shoulder_work_started=False)
    (OUT/'elbow_handoff_provenance.json').write_text(json.dumps(evidence,indent=2),encoding='utf8')
    prior=OUT/'CURRENT_CANDIDATE.json'
    if prior.exists() and not (OUT/'CURRENT_CANDIDATE_M10_R1.json').exists():
        (OUT/'CURRENT_CANDIDATE_M10_R1.json').write_bytes(prior.read_bytes())
    index=dict(version='M10-R1 rigid + continuous-channel TPU trial candidate',
        scope='elbow-only design handoff, no further shoulder design in this request',
        review=str((OUT/'ELBOW_HANDOFF.md').relative_to(ROOT)),
        assembly=str((CONT/'assembly.json').relative_to(ROOT)),
        neutral_MuJoCo=str((CONT/'neutral_view.xml').relative_to(ROOT)),
        neutral_image=str((CONT/'elbow_neutral.png').relative_to(ROOT)),
        rigid_CAD=str((BASE/'cad/paired_swivel_arm_review.step').relative_to(ROOT)),
        rigid_motion=str((BASE/'rigid_motion.gif').relative_to(ROOT)),
        coupons=str((OUT/'elbow_test_coupons/manifest.json').relative_to(ROOT)),
        declared_rigid_path_degrees=[0,150],actual_TPU_driven_range_degrees=None,
        TPU_actuation_verified=False,complete_assembly_verified=False,physical_test_verified=False,
        frozen=False,shoulder_muscle_work_authorized_next=False,
        pending=['printed coupon calibration','full muscle deformation and wire friction','passive pivot equilibrium and installed usable range','crimp, preload, tolerances, retention and complete physical assembly'])
    prior.write_text(json.dumps(index,indent=2),encoding='utf8')
    (OUT/'CURRENT_CANDIDATE.md').write_text('''# 当前肘部设计阶段交付

入口：[肘部设计与验证边界](ELBOW_HANDOFF.md)。本轮在肘部停止，不开展肩部肌肉设计。

M10-R1刚性结构＋两块真实连续通道TPU候选；127件名义装配、125实体刚性STEP，四类打印对比样件。
当前MuJoCo零位图包含实际TPU制造网格；运动动画只展示规定刚性路径。

0–150°不是已实现的肌肉驱动范围；TPU驱动、被动平衡、整机装配和实物性能仍未验证。
交付为设计候选，未制造定版。旧版本及失败诊断保留。
''',encoding='utf8')
    print('elbow-only handoff indexed; interfaces',len(interfaces),'actuation unverified',flush=True)


if __name__=='__main__':run()
