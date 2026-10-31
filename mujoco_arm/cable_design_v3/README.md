# Cable design V3

V3 is a route-repair version derived from V2. V2 files and the interactive demo
remain unchanged.

## Build

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v3\build_v3.py
```

Run the V3 checker:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\cable_design_v2\check_v2.py `
  --model D:\XJTLU\PSP\mujoco_arm\cable_design_v3\arm_cable_design_v3.xml `
  --config D:\XJTLU\PSP\mujoco_arm\cable_design_v3\cable_design_v3.json `
  --matrix D:\XJTLU\PSP\mujoco_arm\cable_design_v3\action_cable_matrix_v3.json `
  --evaluation D:\XJTLU\PSP\mujoco_arm\exports\cable_design_v3\evaluation_v3.json `
  --out D:\XJTLU\PSP\mujoco_arm\exports\cable_design_v3\check_report_v3.json
```

## V3 changes

- K7 GH mandatory cable changed from G01 to repaired G02.
- G02 guide/anchor set remains wrapped at all checked K2/K7 fractions.
- S03 anchors changed for SLIDE_BACK clearance.
- S05/S06 anchors changed for the -5 degree SCAP_INT fallback.

## Results

```text
contraction failures: none
cable clearance failures: none
tension failures: K2, K3, K7, SLIDE_FWD
visual clearance failures: unchanged from V2
```

Important distinction:

- K7 contraction and cable clearance are repaired.
- K7 force closure is not repaired. G02 route motion no longer shrinks too
  much, but the available positive cable tensions cannot balance the K7
  generalized force. The current K7 relative residual is about 0.56 and the
  selected cable set does not positively span the four-dimensional force space.
- Full `-20 deg` SCAP_INT is infeasible with 2 mm cable clearance because the
  S05/S06 scapula-side anchors are themselves inside the thorax. V3 explicitly
  evaluates a `-5 deg` fallback and records the full-target failure.

## Outputs

```text
cable_design_v3.json
action_cable_matrix_v3.json
arm_cable_design_v3.xml
..\exports\cable_design_v3\evaluation_v3.json
..\exports\cable_design_v3\tension_report_v3.json
..\exports\cable_design_v3\route_validation_v3.json
..\exports\cable_design_v3\route_repair_report_v3.json
..\exports\cable_design_v3\check_report_v3.json
..\exports\cable_design_v3\cable_routes_v3.csv
..\exports\cable_design_v3\SolidWorks_hardpoints_v3.csv
```

Design documents:

```text
design/ACTION_CABLE_MATRIX_V3.md
design/CABLE_ANCHORS_V3.md
```
