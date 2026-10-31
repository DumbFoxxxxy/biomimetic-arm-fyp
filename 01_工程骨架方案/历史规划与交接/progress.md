# Progress

## 2026-10-28

### V8 Five-DOF Shoulder Chain

- **Status:** complete; gate blocked
- Read the approved V8 plan and the existing V6/V7 shoulder implementation.
- Confirmed that V7 currently uses only one AC hinge and that its shoulder
  contact gate fails at 42.948 mm detachment and -4.548 mm mesh clearance.
- Selected native AC ball-joint representation and `s,t,phi` surface targets.
- Started planning-file transition from the historical V2 work to V8.
- Implemented `shoulder_chain_v8.json`, surface target generation, native AC
  ball-joint MJCF, direct kinematics and bounded deterministic IK.
- Generated `v8_shoulder_open_chain.xml` and `v8_shoulder_contact.xml`.
- Corrected the zero pose to use the established V6 SC/AC offsets.
- Identified and removed a full-pose overconstraint by defining `phi` as the
  rotation about the R400 outward normal.
- Added permitted surface path corrections for K6, SCAP_EXT, SCAP_INT and the
  axial tests; reversed `s` to the anatomical forward direction.
- Final IK residual: `3.756e-6 mm` position and `1.072e-6 deg` orientation.
- Remaining failure: `SLIDE_BACK_1.00` at `-3.882 mm` full-mesh clearance.
- Best permitted axial correction still gives `-3.876 mm`; maximum feasible
  negative slide under the `-0.5 mm` limit is about `s=-8.5 mm`.
- Contact gate correctly blocked by kinematics; no contact force or frozen
  hardpoint result was manufactured.
- V7 R18 elbow report remains `PASS`; the original 16-cable demo remains
  `PASS`.

#### V8 test results

| Test | Expected | Actual | Status |
| --- | --- | --- | --- |
| Open-chain compile | 3 joints / 5 scalar DOF | 3 joints / 5 scalar DOF | PASS |
| Contact diagnostic compile | Loads and finite | 4 joints / 6 scalar DOF, finite | PASS |
| Surface IK | <= 0.25 mm / 0.25 deg | <= 3.76e-6 mm / 1.08e-6 deg | PASS |
| Full-mesh clearance | >= -0.5 mm | -3.882 mm at `SLIDE_BACK_1.00` | FAIL |
| Contact gate | Runs only after kinematics | Blocked | EXPECTED BLOCK |
| V7 R18 elbow | PASS | PASS | PASS |
| 16-cable demo | PASS | PASS | PASS |

## 2026-10-26

- Started isolated K6/S01 engineering baseline.
- Existing demo files are read-only for this task.
- Generated analytic R400 band and 50 x 130 mm scapula pad.
- Built isolated 4-body K6 model with one active S01 tendon and one guide.
- Implemented search, audit, viewer and pose rendering.
- Full audit currently fails: maximum 30.383 deg at 40 mm pull and 14.850 mm
  full-scapula visual penetration.
- `py_compile` passed for every engineering script.
- Existing 16-cable demo rebuilt and passed its full check:
  `PASS: all demo groups produce finite, visible motion.`
- Added `mujoco_arm/cable_design_v2/` with 16-route model generation, gravity,
  contact band/pad, true guide wrapping, positive-tension QP, action grouping,
  sensitivity, staged reports and SolidWorks hardpoint export.
- V2 structural and minimum-motion checks pass. Full force/clearance checks
  remain failing and are retained in `check_report_v2.json`.
- With group-level cable reassignment, K1-K5, both scapular rotations, the slide
  actions and the pectoral action now have positive-tension solutions. K6/K7
  tension and several clearance checks remain open.
- Final V2 check: structural PASS, minimum-motion PASS, full acceptance FAIL.
  Tension failures are K6/K7; SCAP_EXT exceeds contraction; route and full
  scapula visual clearances remain open.
- Four-muscle optimization replaced S01-S04 anchors with the paper-inspired
  crossed upper/lower layout. K6 and K7 now pass positive spanning and tension.
- Current remaining failures: K7 contraction, K2/SLIDE_BACK/SCAP_INT cable
  clearance, and full scapula visual clearance.
- Added V3 route repair. Contraction and all active cable-clearance checks now
  pass. K2/K3/K7/SLIDE_FWD still have positive-tension regressions, and full
  scapula visual clearance remains unresolved.
- V3 final check: structural PASS, minimum-motion PASS, contraction failures
  none, cable-clearance failures none, force failures K2/K3/K7/SLIDE_FWD.
