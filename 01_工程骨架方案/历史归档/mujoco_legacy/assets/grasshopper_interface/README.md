# Grasshopper interface

Run:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\export_grasshopper_manifest.py
```

The exporter writes:

- `exports/grasshopper_muscles.json`
- `exports/grasshopper_muscles.csv`

Each enabled muscle contains its shape class, cable count, world-space origin,
insertion and via points, and baseline `x/y/z/p` parameters. Grasshopper should
generate one STEP/STL per muscle plus its cable centreline, using the same
muscle id as the filename prefix.

The intended PSM construction order is:

1. Define the muscle envelope from the origin, via and insertion curve.
2. Orient the X-lattice so its contraction direction follows the cable path.
3. Add dense end blocks and wire channels.
4. Add the outer protective shell.
5. Export `psm_<muscle_id>.step` and `psm_<muscle_id>_cable.dxf`.
