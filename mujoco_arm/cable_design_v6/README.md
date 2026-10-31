# V6 continuous-trajectory cable arm

V6 is an isolated continuation of V5. It does not modify V2-V5, the 16-cable
demo, or any SolidWorks file.

## Current status

```text
joints:                  7
cable channels:          17
positive spanning gate:  disabled
static action pass:      false
elbow antagonist pass:   false
trajectory full pass:    false
frozen hardpoints:       false
```

The V6 model implements:

- two independent elbow cables with real R10 wrap;
- a single-sided R400 shoulder roller with compliant preload metadata;
- three ADL trajectories and one all-action trajectory;
- task-conditioned nonnegative tension feasibility;
- 17/18/19/20-channel sweep;
- Chinese trajectory playback with wrist-path and error plots.

The current hardpoints are diagnostic only. Do not use them as final
manufacturing dimensions until `check_report_v6.json` reports
`frozen_hardpoints: true`.

## Commands

```powershell
D:\MujocoEnv\Scripts\python.exe build_v6.py
D:\MujocoEnv\Scripts\python.exe audit_v6.py
D:\MujocoEnv\Scripts\python.exe export_v6_docs.py
D:\MujocoEnv\Scripts\python.exe sweep_channels_v6.py
D:\MujocoEnv\Scripts\python.exe run_v6_gui.py --check
D:\MujocoEnv\Scripts\python.exe run_v6_gui.py
```

## Main outputs

```text
cable_design_v6.json
action_cable_matrix_v6.json
trajectory_library_v6.json
arm_cable_design_v6.xml
../exports/cable_design_v6/check_report_v6.json
../exports/cable_design_v6/trajectory_report_v6.json
../exports/cable_design_v6/channel_sweep_v6.json
../exports/cable_design_v6/SolidWorks_hardpoints_v6.csv
```

