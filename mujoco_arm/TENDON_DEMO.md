# 16-tendon interactive MuJoCo demo

This model is an interactive demonstration of the current 7-DOF arm:

```text
scapula slide + scapula up rotation + scapula external rotation
+ GH flexion + GH abduction + GH axial rotation
+ elbow flexion
```

The 16 cable paths use the final A/G/B coordinates. Each cable is a MuJoCo
spatial tendon and has one position actuator. A lower target length pulls the
cable; the arm moves because of the resulting tendon forces, not because a
joint-angle motor is driving the pose.

The visual meshes are exported from the current SolidWorks parts:

```text
Thorax.SLDPRT
clavicle.SLDPRT
scapula.SLDPRT
humerus.SLDPRT
ulna.SLDPRT
```

The clavicle now has the SC two-axis joint to the thorax. The scapula is
attached to the clavicle through the AC external-rotation hinge. This gives the
intended three scapular DOF while keeping the clavicle motion visible.

## Run the control panel

Double-click:

```text
open_tendon_demo_gui.bat
```

Or run:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_tendon_demo_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_tendon_demo_gui.py
```

The Tk panel is in Chinese, with 16 cable-pull sliders in millimetres,
Chinese cable labels, pose preset buttons and a current-preset status line.
The MuJoCo viewer shows the bones and tendon paths.

After changing the SolidWorks parts, regenerate the meshes with:

```powershell
powershell -ExecutionPolicy Bypass -File D:\XJTLU\PSP\mujoco_arm\export_tendon_demo_meshes.ps1
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\prepare_tendon_demo_meshes.py
```

Preset buttons:

```text
K1 Elbow
K2 GH flex
K3 GH abd
K4 GH internal
K5 GH external
K6 Scap up
K7 Combined
Scap down
Slide +
Slide -
Scap ext +
Scap ext -
胸大肌抱合
```

## Run the viewer without the extra panel

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_tendon_demo.py
```

In this version, use the actuator sliders in the MuJoCo control panel. Slider
values are absolute target tendon lengths in metres. Reduce a value to pull the
wire. Number keys `1`-`9` and letter keys `a`-`c` load presets.

## Validate

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_tendon_demo.py
```

The report is written to:

```text
mujoco_arm/exports/tendon_demo/check_report.json
```

Current representative results:

| Test | Result |
| --- | --- |
| E01 pull 40 mm | elbow 120.0 deg |
| G03 pull 40 mm | GH abduction 44.9 deg |
| S01 pull 40 mm | scapula upward rotation 33.5 deg |
| K7 G01 + S01 | combined GH and scapula motion |

## Render representative poses

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\render_tendon_demo.py
```

Images are written to:

```text
mujoco_arm/exports/tendon_demo/zero.png
mujoco_arm/exports/tendon_demo/K1_E01.png
mujoco_arm/exports/tendon_demo/K3_G03.png
mujoco_arm/exports/tendon_demo/K6_S01.png
mujoco_arm/exports/tendon_demo/P_both.png
```

## Current model assumptions

- Bone visuals use the current SolidWorks meshes.
- The thorax, clavicle and scapula have collision proxies; the scapula is kept
  out of the thorax by contact constraints.
- The SC elevation is limited to `-5..30 deg`, SC protraction to `-15..15 deg`
  and AC external rotation to `-15..15 deg` so K6/K7 cannot drive the
  scapula through the thorax.
- K6/K7 tested contact penetration is below 0.35 mm. If the scapula appears
  to pass through the transparent thorax shell in the viewer, inspect from the
  dorsal camera and use the contact report as the physical check.
- The S02 scapular downward path was retuned to
  `T120_160 -> superior scapular corner` for a clearer downward/medial pull.
- Tendons are ideal piecewise-linear spatial tendons.
- PSM lattice, cable mass, friction and bending stiffness are not simulated.
- Cable actuators are pull-only using a negative force range.
- Gravity is off by default so the arm stays in the zero pose until a cable is
  pulled.
- This is a proof-of-principle demonstrator, not the final physical validation.
