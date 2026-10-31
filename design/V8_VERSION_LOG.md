# V8 version log

## V8.0 - 2026-10-28

- Created an isolated V8 shoulder directory and configuration.
- Replaced the single AC hinge with a native MuJoCo ball joint.
- Kept the original SC two-axis universal joint and 150 mm clavicle.

## V8.1 - 2026-10-28

- Implemented `s,t,phi` surface targets and deterministic five-DOF IK.
- Initial full-pose orientation constraints produced a maximum position error
  of 2.728 mm and a minimum full-mesh clearance of -12.349 mm.
- Diagnosed the full orientation target as over-constraining the five-DOF chain.

## V8.2 - 2026-10-28

- Changed `phi` to constrain rotation about the R400 outward normal while
  leaving the remaining chain redundancy to geometry.
- IK position and orientation residuals fell below `4e-6 mm` and `2e-6 deg`.
- The uncorrected K6 path still reached -13.132 mm full-mesh clearance.

## V8.3 - 2026-10-28

- Reversed the `s` direction to the anatomical forward direction.
- Added permitted surface-path corrections: K6, SCAP_EXT, and AXIAL now combine
  their primary angle/translation with small `s/t` shifts.
- All requested targets passed except `SLIDE_BACK_1.00`.
- Minimum full-mesh clearance improved to `-3.882 mm` at `SLIDE_BACK`.
- Axial correction search did not remove the failure.

## V8.4 - 2026-10-28

- Added the contact-gate prerequisite, regression checks, report export and
  shoulder-gate documentation.
- Contact evaluation is correctly blocked by kinematics.
- V7 R18 elbow regression: PASS.
- 16-cable demo regression: PASS.
- Frozen hardpoints: false.
- Largest clearance-feasible negative slide under the current thorax:
  `s=-8.5 mm` at about `-0.430 mm` mesh clearance.

