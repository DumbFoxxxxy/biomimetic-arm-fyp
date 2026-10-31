# Findings and Decisions

## V8 Requirements

- Preserve SC `(40,235,-106)`, AC `(156.86,233.24,-11.97)` and `SC-AC=150 mm`.
- Use two SC hinges, a rigid clavicle and a native three-DOF AC ball joint.
- Define shoulder motion with `s` along R400, `t` along the cylinder axis and
  `phi` around the local outward normal.
- Pass kinematics and full-mesh clearance before testing any roller contact.
- Do not modify V2-V7, the 16-cable demo or SolidWorks.

## Existing Model Findings

- V6/V7 model `clavicle_elevation -> clavicle_protraction -> clavicle_body` with
  only one AC hinge named `ac_ext`.
- V7 contact failure reached 42.948 mm roller detachment and -4.548 mm full-mesh
  clearance under the existing three-hinge shoulder coordinates.
- V7's elbow R18 gate already passes and must remain unchanged.
- The R400 contact band uses radius 400 mm, x range 55-140 mm and y range
  60-280 mm in design coordinates.
- Existing scapula surface landmarks are P `(62.18,210.10,91.39)`, SUP
  `(72.76,268.45,82.21)` and INF `(44.55,112.87,106.69)`.

## Technical Decisions

| Decision | Rationale |
| --- | --- |
| Native `ball` joint | Exact AC spherical joint and robust quaternion state |
| `s,t,phi` surface frame | Separates desired shoulder motion from five internal joints |
| Two-stage open-chain/contact models | Prevents contact forces from masking an unreachable kinematic target |
| Deterministic IK from previous sample | Avoids branch switching and non-repeatable results |
| Constrain `phi` only, not full pad orientation | `s,t` already determine the contact point; a full frame imposes a sixth independent constraint on five DOFs |
| Reuse V6 neutral offsets as the V8 zero | Prevents an artificial identity-orientation shoulder pose |
| Use clearance-guided surface path corrections | Allows `s/t` motion to compensate for the curved thorax without relaxing penetration limits |

## V8 Results

- Open-chain model: 3 joint objects, 5 scalar DOFs, `nq=6`, `nv=5`.
- Contact diagnostic model: 4 joint objects, `nq=7`, `nv=6`.
- Maximum IK position residual: `3.756e-6 mm`.
- Maximum IK `phi` residual: `1.072e-6 deg`.
- K6, SCAP_EXT, SCAP_INT, SLIDE_FWD, AXIAL and AXIAL_NEGATIVE pass the
  full-mesh clearance gate.
- `SLIDE_BACK_1.00` fails at `-3.882 mm` mesh clearance.
- Best axial correction for the failed endpoint is only `-3.876 mm`.
- Largest clearance-feasible negative slide under `-0.5 mm` is about
  `s=-8.5 mm` with `-0.430 mm` clearance.
- The contact gate is blocked by kinematics, so no roller force result is
  claimed and `frozen_hardpoints=false`.

## Literature Note

- Blache 2024 per-movement results: median 0.38-0.52 mm, P95 1.4-1.5 mm,
  maximum 3.8-4.1 mm.
- The 10.4-11.4 mm maximum applies when one ellipsoid is shared across
  different movements.
- The 1-5 mm allowance is a project engineering tolerance, not a paper limit.

## Errors and Resolutions

| Issue | Resolution |
| --- | --- |
| Full orientation target made five-DOF IK infeasible | Reduced the target to contact point plus `phi` projection |
| Zero target contradicted the established shoulder neutral | Used V6 zero-pose offsets |
| One slide direction penetrated the R400 envelope | Searched allowed axial corrections and reported the remaining geometric limit |

