# Cable-drive exports

The cable-drive model needs two inputs only:

1. Local binary STL files for `thorax`, `clavicle`, `scapula`, `humerus`
   and `ulna`.
2. Named SolidWorks coordinate systems `SC_J`, `ACST_J`, `GH_J` and
   `ELBOW_J`.

No STEP, OBJ, muscle mesh or cable mesh is required for this version.

## Two ways to describe the pose

The recommended way is a coordinate system on *both* bones of every joint:

```text
SC_J     clavicle   + thorax
ACST_J   scapula    + clavicle
GH_J     humerus    + scapula
ELBOW_J  ulna       + humerus
```

Both frames of a pair must sit at the same anatomical joint centre with the
same axis meanings. The exporter then places each child bone so its frame
coincides with the parent frame, and the pose of the parts inside
`biomimetic_structure.SLDASM` stops mattering.

If the parent-side frame is missing, the exporter falls back to the pose the
part happens to have in `biomimetic_structure.SLDASM`.

Check the export before building the model:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_bone_layout.py
```

## SolidWorks preparation

Create the four coordinate systems at the joint centres, either at assembly
level or inside the owning bone part. Part-level frames are recommended when
the assembly components are not in a clean neutral pose. Their axes are
interpreted in this order:

```text
SC_J:    X abduction/adduction, Y elevation/depression, Z axial rotation
ACST_J:  X upward rotation, Y posterior tilt, Z external rotation
GH_J:    X rotation, Y flexion, Z abduction
ELBOW_J: X flexion/extension
```

Recommended owning parts:

```text
SC_J:    clavicle
ACST_J:  scapula
GH_J:    humerus
ELBOW_J: ulna
```

Run the read-only export macro:

`D:\XJTLU\PSP\model\validation\ExportCableDriveData.swb`

It writes:

- `bones_local/*.stl`
- `component_transforms.txt`
- `joint_frames.txt`
- `part_frames.txt`

## Python preparation

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\prepare_cable_drive_exports.py
```

This applies the assembly transforms, converts the meshes to metres and writes:

- `bones/*.stl`
- `joint_frames.json`

The MuJoCo builder stops with an explicit error if any of these files are
missing or if a joint frame is still marked `MISSING`.
