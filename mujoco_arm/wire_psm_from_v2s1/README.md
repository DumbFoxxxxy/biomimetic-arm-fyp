# V2-S1 钢丝与TPU肌肉接口

当前入口为 `../exports/wire_psm_from_v2s1/PSM-W1/README.md`（评审候选，不是加工定版）。原骨架资产不修改。

选择的结构是A固定肌腹端、C移动肌腹端、G偏转滑轮固定进线切点、B另一骨骼止点。电机主动钢丝穿过TPU至C终止；另一个固定长度的肌腱从C经滑轮到B。不要把C当骨孔，不要把G当实际轮心。

当前流程：`design_modules.py` → `solve_swivel_interfaces.py` → `audit_interfaces.py` → `export_layout.py` → `write_report.py`。`render_layout.py`用带Pillow的Python运行。`export_native_layout.ps1`通过已验证的强类型SolidWorks桥创建独立参考草图零件，不改源CAD。

`audit_seed.py`为V2原走线复查；`fit_mounts.py`为早期仅法向偏置试验；`solve_belly_interfaces.py`为已否决的自由球形导向耦合试验。它们的旧输出不能作为当前偏转滑轮接口的证明。实际钢丝圆弧使用common.swivel_route，XML中的球形标记仅用来保存G位置。

仅评估规定姿态几何、肌腹端行程及明确列出的空间/应变情景。未完成支架实体、外部护套线束、肌肉间形变相容性、实际TPU能力或电机驱动验收。
