# MuJoCo reference arm

This folder contains the first kinematic model for the upper arm.

For the new five-bone, cable-and-slider shoulder/elbow simulation, see
`CABLE_DRIVE.md`. It uses the new `model/exports` data and does not reuse the
old radius/pronation model.

## Run

Use the existing virtual environment:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_pronation.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_scapula_rhythm.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\validate_contract.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_arm_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_psm_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_muscle_paths.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\export_grasshopper_manifest.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_viewer.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_shoulder_rhythm.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_psm_viewer.py
```

## What this first model contains

- BodyParts3D reference meshes for the clavicle, scapula, humerus, radius and ulna.
- Glenohumeral motion represented by three independent hinge joints:
  `shoulder_rotation`, `shoulder_flexion` and `shoulder_abduction`.
- The ulna carries the elbow hinge and remains the forearm reference bone.
- The proximal radiocapitellar contact point is the pivot and remains fixed
  while the radius pronates around the line from that contact point to the
  distal ulnar head.
- The sternoclavicular joint has elevation/depression, protraction/retraction
  and axial rotation hinges.
- The acromioclavicular/scapulothoracic motion has upward rotation, posterior
  tilt and external rotation hinges.
- `run_shoulder_rhythm.py` couples those hinges to shoulder flexion and
  abduction with a first-pass scapulohumeral rhythm.
- Eleven position actuators so every scalar joint can be controlled from the
  GUI sliders.
- Gravity disabled for the first kinematic pass.

The config-driven musculoskeletal model is built by:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\validate_contract.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_arm_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_psm_model.py
```

It uses:

- `config/arm_model.json` for bones, joints, thoracic references and ligaments.
- `config/muscle_manifest.json` for the muscle candidate pool and PSM parameters.
- `arm_psm_generated.xml` as generated output.

The generated model currently enables the seven required muscles:
supraspinatus, infraspinatus, anterior/middle/posterior deltoid, elbow
flexors, and triceps. The remaining muscles are present in the manifest but
disabled until optimization enables them.

The optimizer scaffold is:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\optimize_muscles.py --trials 20
```

It uses Optuna NSGA-II when available and a deterministic random fallback
otherwise. Results are written to `exports/pareto_designs.json`.

## Actual CAD exports

The SolidWorks exports are adopted with:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\import_cad_export.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_arm_model.py `
  --arm-config D:\XJTLU\PSP\mujoco_arm\config\arm_model_cad.json `
  --muscle-config D:\XJTLU\PSP\mujoco_arm\config\muscle_manifest_cad.json
```

Generated CAD-based files are:

- `arm_cad_reference.xml`
- `arm_psm_generated_cad.xml`
- `config/arm_model_cad.json`
- `config/muscle_manifest_cad.json`
- `config/cad_registration.json`

Run the CAD model viewer with:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_psm_viewer.py `
  --model D:\XJTLU\PSP\mujoco_arm\arm_psm_generated_cad.xml
```

## Rhino PSM generation

Open Rhino 8 and run:

```text
_-RunPythonScript "D:\XJTLU\PSP\mujoco_arm\grasshopper_interface\generate_psm_rhino.py"
```

The script reads `exports/grasshopper_muscles_cad.json` and exports:

- `exports/rhino_psm/right_arm_psm.3dm`
- `exports/rhino_psm/right_arm_psm.stl`

The CAD registration is generated from three anatomical landmarks per bone.
The current sample residuals are approximately 0.6-12.5 mm.

The clavicle and scapula joints are simplified hinge representations of the
real sternoclavicular and acromioclavicular joints. Replace their axes and
ranges after measuring the real joint surfaces and ligaments.

## Values that must be replaced from your own assembly

The current coordinates are derived from the BodyParts3D meshes in
`D:\XJTLU\PSP\骨骼`. Replace these pivot values with measurements from
`part.SLDASM`:

- sternoclavicular joint
- acromioclavicular joint
- glenohumeral joint centre
- elbow centre and elbow axis
- wrist centre and forearm rotation axis

Also replace the current joint ranges:

- shoulder rotation: `-90` to `90` degrees
- shoulder flexion: `-45` to `180` degrees
- shoulder abduction: `-45` to `180` degrees
- elbow flexion: `0` to `145` degrees
- forearm pronation: `-90` to `90` degrees

These are placeholder limits for checking the model, not final anatomical
limits.

## Next replacement step

Export the left-arm parts from SolidWorks as ASCII OBJ or STL with their
assembly transforms applied, replace the files in `assets`, and update the
four joint centres. The XML structure can otherwise remain unchanged.
