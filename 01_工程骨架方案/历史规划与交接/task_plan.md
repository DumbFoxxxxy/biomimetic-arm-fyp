# Task Plan: V8 Five-DOF Shoulder Chain

## Goal

Build an isolated five-DOF scapular mechanism with a two-axis SC universal
joint, a rigid clavicle and a native AC ball joint, then gate the mechanism on
surface-frame inverse kinematics and full-mesh clearance before adding R400
roller contact.

## Current Phase

Complete - implementation finished with the kinematic gate blocked

## Phases

### Phase 1: Requirements and Existing Model Audit
- [x] Read the approved V8 plan
- [x] Audit the V6/V7 shoulder tree, configurations and contact code
- [x] Confirm V8 is isolated from V2-V7 and the demo model
- **Status:** complete

### Phase 2: V8 Kinematic Model and Surface Targets
- [x] Create the V8 configuration, target generator and MJCF builder
- [x] Implement the native AC ball joint and five-DOF IK
- [x] Generate open-chain and contact XML entrypoints
- **Status:** complete

### Phase 3: Kinematic and Mesh-Clearance Gate
- [x] Implement deterministic IK continuity and residual checks
- [x] Implement analytic scapulothoracic clearance checks
- [x] Run the kinematic gate and emit `kinematics_report_v8.json`
- **Status:** complete, gate failed on `SLIDE_BACK_1.00`

### Phase 4: R400 Contact Gate
- [x] Add single and dual passive roller candidates
- [x] Add prerequisite blocking by the kinematic gate
- [x] Emit `contact_report_v8.json`
- **Status:** complete, correctly blocked by kinematics

### Phase 5: Integration, Regression and Documentation
- [x] Apply the overall gate and emit `integration_report_v8.json`
- [x] Run V7 elbow and 16-cable demo regression checks
- [x] Write `V8_SHOULDER_CHAIN.md` and final progress records
- **Status:** complete

## Acceptance

- SC-AC center distance remains 150 mm within 0.05 mm.
- IK position residual is <= 0.25 mm and orientation residual <= 0.25 deg.
- Full scapula mesh clearance is >= -0.5 mm.
- Contact roller detachment is <= 1 mm and penetration is <= 0.5 mm.
- Nominal roller contact force is 3-10 N and preload sensitivity is stable.
- No cable optimization or frozen hardpoint export occurs before both gates pass.

## Decisions

| Decision | Rationale |
| --- | --- |
| Native MuJoCo ball joint at AC | Exact three-DOF spherical motion without serial-axis singularities |
| Surface targets `s,t,phi` | Keeps net shoulder behavior on the R400 engineering surface |
| K6 and SCAP_EXT/INT merged into `phi` | Required by the approved strict three-coordinate output |
| Kinematic gate before contact | V7 showed contact tuning can hide a closed-chain geometry error |
| Leave non-phi tilt to the mechanism redundancy | A full pose target over-constrained the five-DOF chain |

## Errors

| Error | Attempt | Resolution |
| --- | --- | --- |
| Full orientation target produced 2.728 mm residual | 1 | Constrained only the R400 normal rotation `phi` |
| Initial zero pose forced world identity and caused 18.039 mm penetration | 1 | Reused the established V6 neutral SC/AC offsets |
| Positive `s` initially moved opposite the anatomical forward direction | 1 | Reversed the R400 arc direction |
| `SLIDE_BACK=-15 mm` penetrates -3.882 mm | 1 | Kept the threshold; reported maximum feasible `-8.5 mm` and blocked contact |

