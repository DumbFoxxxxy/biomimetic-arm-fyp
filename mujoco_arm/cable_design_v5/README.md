# V5 cable-design baseline

V5 is an isolated branch of the cable-design work. It does not modify V2, V3,
V4, `arm_tendon_demo.xml` or any SolidWorks file.

## Current status

```text
joints:             7
channels tested:    16, 17, 18 (G08 and S07)
current baseline:   17
full pass:          false
frozen hardpoints:  false
```

The full scapula mesh now passes the 21-pose clearance check. The remaining
blocking issue is mechanical force feasibility: several actions cannot satisfy
the hard wrench-closure/positive-spanning condition, and the current elbow
route cannot produce opposing E01/E02 moments over the complete stroke.

Do not use the V5 hardpoints as frozen manufacturing data until
`check_report_v5.json` reports `frozen_hardpoints: true`.

## Main commands

```powershell
D:\MujocoEnv\Scripts\python.exe build_v5.py
D:\MujocoEnv\Scripts\python.exe audit_v5.py
D:\MujocoEnv\Scripts\python.exe export_v5_docs.py
```

Optimization and diagnosis tools:

```text
optimize_pose_v5.py
optimize_e_routes_v5.py
optimize_s_routes_v5.py
sweep_channels_v5.py
diagnose_e_guides_v5.py
```

## Main outputs

```text
arm_cable_design_v5.xml
cable_design_v5.json
action_cable_matrix_v5.json
../exports/cable_design_v5/check_report_v5.json
../exports/cable_design_v5/channel_sweep_v5.json
../exports/cable_design_v5/capstan_report_v5.json
../exports/cable_design_v5/contact_audit_v5.json
../exports/cable_design_v5/ellipsoid_audit_v5.json
../exports/cable_design_v5/SolidWorks_hardpoints_v5.csv
../../design/CABLE_ANCHORS_V5.md
../../design/ACTION_CABLE_MATRIX_V5.md
../../design/V5_DESIGN_DECISIONS.md
```
