# K6/S01 engineering baseline

This directory is isolated from the interactive 16-cable demo. It contains one
active S01 spatial tendon, the existing three shoulder hinges, analytic contact
geometry, and a thorax-mounted cylindrical guide.

## Build and audit

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\engineering\generate_contact_assets.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\engineering\build_k6_baseline.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\engineering\k6_audit.py
```

The audit exits non-zero when the baseline fails. The machine-readable result is:

```text
mujoco_arm/exports/k6_baseline/report.json
```

## Inspect and render

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\engineering\run_k6_baseline.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\engineering\render_k6_baseline.py
```

Rendered poses are written to:

```text
mujoco_arm/exports/k6_baseline/poses/
```

## Current result

The selected candidate is the best true-wrap route found inside the agreed
search bounds. It does not pass the baseline:

| Metric | Result |
| --- | ---: |
| Zero-pose contact gap | 0.348 mm |
| Maximum upward rotation at 40 mm pull | 30.383 deg |
| Maximum pad penetration | 0 mm |
| Maximum free separation | 26.303 mm |
| Maximum full-scapula visual penetration | 14.850 mm |
| Minimum S01 path clearance to thorax | 6.062 mm |

The result means the current three-hinge mechanism plus one thorax-mounted
guide cannot meet the 36 deg requirement within the specified search bounds,
and the existing scapula visual mesh penetrates the analytic thorax envelope.
The search bounds were not expanded after this failure.

## Search records

```text
assets/route_search.json
assets/route_search_active.json
```

The random active-wrap search is the more relevant record because it rejects
guides that do not actually enter the cable wrap.
