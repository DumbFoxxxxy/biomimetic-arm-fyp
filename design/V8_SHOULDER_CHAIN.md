# V8 five-DOF shoulder chain

## Purpose

V8 replaces the single AC hinge with an exact three-DOF AC ball joint
and keeps the SC two-axis universal joint and the 150 mm clavicle.
Shoulder motion is defined by the surface coordinates `s`, `t`, and
`phi`; cable optimization is outside this gate.

## Model

- SC center: `(40, 235, -106) mm`
- AC center: `(156.86, 233.24, -11.97) mm`
- SC-AC design length: `150 mm`
- MuJoCo joints: `sc_elevation`, `sc_protraction`, `ac_ball`
- Scalar shoulder DOFs: `5`
- Full scapula and clavicle meshes remain visual-only

## Gate Results

| Check | Result |
| --- | --- |
| Kinematic gate | `False` |
| Contact gate | `False` |
| Integration gate | `False` |
| Frozen hardpoints | `False` |
| Maximum IK position residual | `3.75636e-06 mm` |
| Maximum IK orientation residual | `1.07222e-06 deg` |
| Minimum full-mesh clearance | `-3.882 mm` |
| Minimum contact-pad clearance | `0.285 mm` |
| Worst target | `SLIDE_BACK_1.00` |

The contact gate reports:

```text
blocked_by_kinematics, kinematics:full_mesh_clearance
```

The maximum currently clearance-feasible negative slide is
`s=-8.5 mm`, clearance `-0.430 mm`.

## Interpretation

The five-DOF chain can follow the new surface targets to micrometre
level. The remaining failure is geometric: at the requested
`SLIDE_BACK=-15 mm`, the rigid scapula mesh intersects the current
R400/ellipsoid thorax envelope even after the permitted surface path
adjustment. Per the approved gate, the R8 roller contact stage was not
used to hide this error and no hardpoint table was frozen.

## Outputs

- `mujoco_arm/cable_design_v8/v8_shoulder_open_chain.xml`
- `mujoco_arm/cable_design_v8/v8_shoulder_contact.xml`
- `mujoco_arm/exports/cable_design_v8/kinematics_report_v8.json`
- `mujoco_arm/exports/cable_design_v8/contact_report_v8.json`
- `mujoco_arm/exports/cable_design_v8/integration_report_v8.json`

The V7 R18 elbow regression and the 16-cable demo regression both pass.
V5/V6/V7 source and model files were not edited.
