# V8 five-DOF shoulder chain

V8 is isolated from V2-V7, the 16-cable demo and SolidWorks.

## Gate order

```text
five-DOF kinematics
  -> surface-frame IK and full-mesh clearance
  -> R8 roller contact
  -> integration
```

The contact gate is not evaluated when the kinematic gate fails.

## Commands

```powershell
D:\MujocoEnv\Scripts\python.exe build_v8_shoulder.py
D:\MujocoEnv\Scripts\python.exe audit_v8_kinematics.py
D:\MujocoEnv\Scripts\python.exe audit_v8_contact.py
D:\MujocoEnv\Scripts\python.exe integration_v8.py
D:\MujocoEnv\Scripts\python.exe export_v8_docs.py
```

Or run the complete gate:

```powershell
.\run_v8_audit.ps1
```

## Current result

The five-DOF model compiles with three shoulder joint objects and five
scalar DOFs. IK reaches the surface targets at micrometre-level residual.

The kinematics gate remains blocked by `SLIDE_BACK_1.00`, whose complete
scapula mesh penetrates the analytic thorax by `3.882 mm`. The contact gate
is therefore intentionally not evaluated.

