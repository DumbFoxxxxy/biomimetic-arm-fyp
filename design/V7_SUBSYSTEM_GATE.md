# V7 subsystem gate decisions

## Elbow sector spool

- Gate: `PASS`
- Baseline radius: `18 mm`
- Minimum effective moment arm: `15.204 mm`
- Maximum active-stroke contraction: `10.562%`
- Opposite moments: `True`
- Real guide wrap: `True`

## Shoulder R400 roller

- Gate: `FAIL`
- Selected roller count: `1`
- Maximum roller detachment: `42.948 mm`
- Minimum visual clearance: `-4.548 mm`
- Maximum roller contact force: `4.298 N`
- Preload sensitivity: `False`

## Integration

- Status: `blocked`
- Gate pass: `False`
- Frozen hardpoints: `False`

## Blocking reasons

- `shoulder:full_mesh_penetration`
- `shoulder:roller_detachment`
- `shoulder:preload_sensitivity`

No V7 SolidWorks hardpoint file was written because the integration gate is closed.
