# Cable design V2

This directory contains the force-aware cable design pipeline. It preserves the
existing 7-DOF mechanism and exports route and hardpoint data for SolidWorks.

## Pipeline

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\prepare_v2.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\build_v2_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\evaluate_v2.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\stage_evaluate.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\sensitivity_v2.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\friction_sensitivity_v2.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\check_v2.py
```

Open the model:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\run_v2.py
```

## Outputs

```text
cable_design_v2.json
action_cable_matrix.json
arm_cable_design_v2.xml
..\exports\cable_design_v2\evaluation_v2.json
..\exports\cable_design_v2\tension_report.json
..\exports\cable_design_v2\route_validation_v2.json
..\exports\cable_design_v2\staged_report.json
..\exports\cable_design_v2\sensitivity_v2.json
..\exports\cable_design_v2\friction_sensitivity_v2.json
..\exports\cable_design_v2\check_report_v2.json
..\exports\cable_design_v2\cable_routes_v2.csv
..\exports\cable_design_v2\SolidWorks_hardpoints_v2.csv
```

## Key finding

The old G coordinates are path bend points, not pulley centres. V2 computes a
true R8 guide centre from the A-G-B geometry. If the bend creates less than 10
degrees of wrap, the point becomes a `pass_through` fairlead rather than a fake
guide.

The S group follows Sodeyama et al.'s four-muscle bladebone layout:

```text
S01/S02: upper scapular points pulled downward
S03/S04: lower scapular points pulled upward from the chest
S05/S06: auxiliary rotation, lift and antagonistic stabilization
```

K6/K7 force S01-S04 to carry at least 2 N each. S05/S06 are optional. A true
positive-spanning LP test is applied to the resulting generalized-force columns.

## Current result

- The model compiles with 7 joints, 16 tendons and real cylindrical wrapping.
- All actions reach their 80 percent pose metric.
- K6 and K7 now have positive-spanning, positive-tension force solutions.
- K6 uses the four mandatory scapular muscles plus S06; its maximum tension is
  about 14 N and contraction is about 17.1 percent.
- K7 adds G01/G03 plus the four mandatory scapular muscles; force closure is
  exact, but G01 contraction exceeds 18 percent and maximum tension is about
  141 N.
- K2, SLIDE_BACK and SCAP_INT still have cable-thorax clearance failures.
- The full scapula visual mesh still penetrates the thorax at some poses.

This is an implemented and exportable design pipeline, but not a frozen physical
layout. Failures remain visible in the machine-readable reports.

Guide friction is approximated by a whole-cable capstan efficiency
`exp(-mu * total_wrap)`. Contact friction is not included in the static
equilibrium QP; `friction_sensitivity_v2.json` records gravity-only settling at
the requested friction values.
