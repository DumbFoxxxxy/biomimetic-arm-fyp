# V6 design decisions

## Architecture

- Seven controlled DOFs, 17-20 independent cable channels.
- Global wrench closure is diagnostic only; local nonnegative tension is checked.
- Two active elbow antagonist cables with R10 guides.
- R400 single-sided crowned POM roller on the scapula with compliant preload.
- Three ADL trajectories plus one all-action sequence.

## Current audit

- Structural pass: `True`
- Static actions pass: `False`
- Elbow antagonism pass: `False`
- Full trajectory pass: `False`
- Frozen hardpoints: `False`

## Trajectory library

- `FRONT_REACH` 前伸取物: 1.6 s
- `HAND_TO_MOUTH` 手到口鼻: 1.4 s
- `LATERAL_RAISE` 侧向上举: 1.8 s
- `ALL_ACTIONS` 全部动作连续演示: 20.0 s

## Failure reasons

- `static:K1`
- `static:K2`
- `static:K3`
- `static:K4`
- `static:K5`
- `static:K6`
- `static:K7`
- `static:SLIDE_BACK`
- `static:SCAP_EXT`
- `static:PECTORAL`
- `elbow:moment_arm`
- `elbow:contraction`
- `trajectory:FRONT_REACH`
- `trajectory:HAND_TO_MOUTH`
- `trajectory:LATERAL_RAISE`
- `trajectory:ALL_ACTIONS`

## Shoulder preload

- Nominal: 5.0 N
- Stiffness: 2.0 N/mm
- Travel: 4.0 mm
