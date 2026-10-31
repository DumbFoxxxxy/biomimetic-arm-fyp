# V7 isolated subsystem validation

V7 does not modify V2-V6, the 16-cable demo, or SolidWorks.

## Gate order

```text
elbow sector spool
  -> shoulder R400 passive roller
  -> K1/K2 integration
```

No integrated arm or SolidWorks hardpoint table is produced unless both local
mechanism stages pass.

## Current result

```text
elbow R18:          PASS
shoulder 1 roller:  FAIL
shoulder 2 rollers: FAIL
integration:        BLOCKED
frozen hardpoints:  false
```

The shoulder result is not caused by the number of rollers. The existing
SC/AC motion and R400 surface separate by more than the allowed 4 mm passive
travel in several prescribed poses.

## Commands

```powershell
D:\MujocoEnv\Scripts\python.exe build_elbow_standalone.py
D:\MujocoEnv\Scripts\python.exe audit_elbow_v7.py
D:\MujocoEnv\Scripts\python.exe build_shoulder_standalone.py
D:\MujocoEnv\Scripts\python.exe audit_shoulder_v7.py
D:\MujocoEnv\Scripts\python.exe integration_v7.py
D:\MujocoEnv\Scripts\python.exe export_v7_docs.py
```

## Main outputs

```text
v7_elbow_standalone.xml
v7_shoulder_standalone.xml
../exports/cable_design_v7/elbow_report_v7.json
../exports/cable_design_v7/shoulder_report_v7.json
../exports/cable_design_v7/integration_report_v7.json
../../design/V7_SUBSYSTEM_GATE.md
```
