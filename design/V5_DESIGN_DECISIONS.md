# V5 design decisions and audit

## Scope

- Seven controlled DOFs and 16-18 independent cable channels.
- K7 target is 126 deg; the previous 150 deg target is retained only as legacy metadata.
- SCAP_INT is formally -5 deg; the full -20 deg target is recorded as infeasible.
- R400 is the manufacturing contact surface; action-specific ellipsoids are sensitivity only.

## Implemented

- Shoulder zero offsets are limited to +/-10 deg and contact pad offset to 0-10 mm.
- K6/K7/SLIDE_BACK/SCAP_INT use five-knot shoulder-path corrections.
- Full scapula vertices are checked at 21 poses against the thorax.
- Tension solving uses exact nonnegative equilibrium, mandatory minimum tension, real guide wrap and hard positive-spanning checks.
- A 16/17/18-channel sweep was generated without overwriting V2/V3/V4.

## Current result

- Structural pass: `True`
- Full pass: `False`
- Frozen hardpoints: `False`

### Channel sweep

| Variant | Channels | Full pass |
| --- | ---: | --- |
| channels16_no_g07 | 16 | False |
| channels17 | 17 | False |
| channels18_g08 | 18 | False |
| channels18_s07 | 18 | False |

## Action status

Passing the complete current V5 audit: K4, SLIDE_FWD, SCAP_INT

Blocking actions: K1, K2, K3, K5, K6, K7, SLIDE_BACK, SCAP_EXT, PECTORAL

K7 achieved: `126.207 deg`; SCAP_INT achieved: `5.002 deg`.

## Blocking findings

- The current E01/E02 elbow topology cannot produce opposing moments over the full 116 deg stroke.
- The preferred and global cable sets still fail hard positive spanning in several actions.
- Adding only G08 or S07 does not remove the K1 failure and does not make the global design pass.
- Therefore this V5 package must not be used as a frozen manufacturing drawing.
