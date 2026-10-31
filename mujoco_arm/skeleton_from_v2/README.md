# V2 骨架优先：实施与交付

唯一选定的几何基线：`../exports/skeleton_from_v2/V2-S1/`。阅读其中 README、acceptance.csv 和 freeze_manifest.json。其他同级导出文件是中间候选或失败试验，不能混用为定版结果。

本阶段保留 V2 七个转动坐标和全部原关节中心，进行实体让位、局部修形和规定运动路径验证。未设计钢丝+TPU驱动，也没有作受力验收。可编辑CAD交付为分面BRep STEP，原生SolidWorks特征历史保留在原目录但未回写。

## 可重现流程

项目根目录运行，几何/仿真使用 `D:/MujocoEnv/Scripts/python.exe`；STEP导出使用Python3.10及项目内CAD依赖；动画使用带Pillow的 `C:/Python314/python.exe`。

1. `build_candidate.py`，`audit_candidate.py`：零位局部修正和初始诊断；它不是最终方案。
2. `swept_reliefs.py --joint-only`：肩/肘局部运动避让。不要用无参数的全胸廓切除试验作为候选，它因破坏肩胛形态已被否决。
3. `fit_clear_shoulder.py --conservative`：寻找保留肩胛主体的肩带路径，输出独立候选路径。
4. `repair_sc_mount.py --fit shoulder_path_fit_conservative.json`，`refine_reverse_clearance.py`：保护轴承与关节区域的支座/反向轮廓修正。
5. `verify_joint_candidate.py --cache reverse_refined_local_meshes.npz --fit shoulder_path_fit_conservative.json`：全部12条路径、首次碰撞边界及回程检查。
6. `export_review.py`，`export_cad_input.py`；Python3.10运行 `export_step.py`、`remesh_step.py`；仿真Python运行 `check_cad_vertices.py`。
7. 动画Python运行 `render_motion.py` 并检查画面；仿真Python运行 `finalize_baseline.py`，证据不齐或哈希不符时会拒绝生成定版包。

脚本路径均相对本目录，运行时提供完整脚本路径。源网格单位为米；几何布尔计算和CAD采用毫米；最终MuJoCo资产为各骨骼本体坐标下的米，STEP为任务书零位全局坐标下的毫米。

## 已否决/仅诊断工具

`inspect_repairs.py`检查让位规律；`repair_sc_corners.py`为未采用的叉耳切角试验；`cad_dispatch.ps1`为只读SolidWorks接口诊断。它们不构成当前交付流程。
