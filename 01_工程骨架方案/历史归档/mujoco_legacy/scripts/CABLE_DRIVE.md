# Cable-driven shoulder and elbow model

This is the first cable-only simulation stage:

- five bones: `thorax`, `clavicle`, `scapula`, `humerus`, `ulna`
- no radius, wrist, hand or forearm pronation
- 10 skeleton DOF plus 8 world-fixed slider DOF
- eight spatial tendons controlled by slider displacement
- fixed, bone-mounted wire guide fixtures
- temporary scapulohumeral rhythm on the clavicle and scapula

## Required SolidWorks data

Create these coordinate systems either at assembly level or inside the owning
part. Part-level coordinate systems are recommended when the assembly parts
are not currently in a clean neutral pose:

```text
SC_J    X abduction/adduction, Y elevation/depression, Z axial rotation
ACST_J  X upward rotation, Y posterior tilt, Z external rotation
GH_J    X rotation, Y flexion, Z abduction
ELBOW_J X flexion/extension
```

Recommended owning parts:

```text
SC_J     clavicle
ACST_J   scapula
GH_J     humerus
ELBOW_J  ulna
```

The export macro first looks for assembly-level frames, then for frames inside
each part. Part frames are combined with the current component transform
automatically.

### Where the bone pose comes from

`prepare_cable_drive_exports.py` supports two sources and picks the better
one automatically.

**Paired frames (preferred).** Give every joint a coordinate system on the
child bone *and* a matching one on the parent bone:

```text
SC_J     clavicle   + thorax
ACST_J   scapula    + clavicle
GH_J     humerus    + scapula
ELBOW_J  ulna       + humerus
```

Each child bone is then placed so that its frame coincides with the parent
frame, which makes the pose inside `biomimetic_structure.SLDASM` irrelevant.
The joint zero position is defined by the frames, so the parts may sit
anywhere in the assembly.

Frame names may be anything that starts with the joint name as long as it is
unambiguous on that bone: `SC_J`, `SC_J_thorax`, `SC_J_P` all resolve.

**Assembly pose (fallback).** If the parent frame is missing, the bone is
placed at the pose recorded in `component_transforms.txt`.

Without the parent-side frames the viewer shows the assembly exactly as it
is, so a scattered assembly gives a scattered arm.

Check the export before building:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_bone_layout.py
```

It prints `PASS` or `FAIL` for two things:

```text
joint frame placement   each frame must sit on the bone that owns it
bone-to-bone gap        bones sharing a joint must be within 15 mm
```

Run it again after every SolidWorks export; it writes
`exports/cable_drive/bone_layout_check.json`.

### The STL export carries its own offset

`ExportCableDriveData.swb` also writes `model/exports/part_probe.txt` with the
volume and centre of mass of every part. The SolidWorks STL exporter writes
the mesh in its own output coordinate system, which is not always the part
origin, so the exported mesh can be shifted by tens of millimetres without
anything in the part looking wrong.

`prepare_cable_drive_exports.py` compares the mesh volume with the
SolidWorks volume and shifts the mesh back onto the part origin using the
centre-of-mass difference. The correction is printed for every bone. If
`part_probe.txt` is missing the correction is skipped, and the bones will
land in the wrong place whenever the STL output coordinate system is not the
part origin.

### Wire anchors live on the bone

`derive_cable_anchors.py` puts every insertion point back on the surface of
its own bone, inside a distance band around the joint it drives, and picks
the side that produces the required moment-arm sign. It writes
`insertion_part_m`, a point in the bone's part coordinates, so the anchor
follows the bone when the assembly is re-posed.

`optimize_guide_points.py` then searches the pulley or guide point that keeps
the moment arm positive over the whole working range.

### Tendon rest length

The tendons use `spring_rest_mode: "reference"`. The rest length is the
neutral-pose length plus `spring_reference_slack_m`, so every cable is slack
with no command and the arm sits still instead of drifting on the balance of
eight taut springs.

MuJoCo tendon springs are two-sided, so a slack cable would push instead of
releasing. The builder therefore adds one `assist_<cable>` general actuator
per channel and the control loop writes `max(0, -k (L - rest))` into it,
which cancels the compressive half and leaves a wire that only pulls.
`cable_drive_common.apply_one_sided_tension` does that; `run_cable_drive.py`
and `test_cable_motions.py` both call it every step.

### Sliders can retract

Each slider travels from `-release_m` to `travel_m`. Positive is a pull,
negative lets the wire go slack. Releasing the idle cables matters: a cable
that stays rigged while another one pulls acts as a spring holding the arm
back. With every idle cable released and only the target channel pulled:

```text
elbow flexion      145.1 deg   target  90
shoulder flexion   191.4 deg   target  90
shoulder abduction 116.3 deg   target  90
shoulder rotation  104.6 deg   target  45
```

Run `test_cable_motions.py` to reproduce; it releases every channel that is
not selected for the motion under test.

Run:

`D:\XJTLU\PSP\model\validation\ExportCableDriveData.swb`

It writes:

```text
model/exports/bones_local/*.stl
model/exports/component_transforms.txt
model/exports/joint_frames.txt
model/exports/part_frames.txt
```

## Convert exports

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\prepare_cable_drive_exports.py
```

This writes metre-based global meshes and `model/exports/joint_frames.json`.

## Generate the cable manifest

The repository already contains `mujoco_arm/config/cable_drive_manifest.json`.
It can be regenerated from the enabled muscle pool:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\make_cable_drive_manifest.py
```

Each muscle uses `guide_mode = constant_moment_arm`. The guide is a rigid
bracket fixture attached to a bone, so it does not have to lie on the bone
surface. Replace these guides with the final muscle paths when the muscle
geometry is frozen. Explicit paths can be supplied with a `guides` list:

```json
{
  "guides": [
    {"body": "scapula_external_rotation", "local_m": [0.0, 0.0, 0.0], "radius_m": 0.006}
  ]
}
```

## Build and inspect

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_cable_drive_model.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\check_cable_drive.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\run_cable_drive.py
```

In the viewer, the `pull_<cable>` actuator sliders are the cable pull inputs.
Positive values are metres of slider extension.

## One-click launch

Two batch files, both double-clickable:

`open_arm.bat` opens the viewer on the model that is already built, so it
starts immediately.  If `arm_cable_drive.xml` is missing it falls back to the
full pipeline once.

`run_cable_drive.bat` is the full pipeline: export the SolidWorks data, measure
the wire anchors, rebuild `arm_cable_drive.xml`, re-optimise the guides, and
open the viewer.  Run this one after the bones change in SolidWorks; it takes
about 40 s because the guide search samples every bone surface.

Use the `pull_<cable>` actuator sliders in the right-hand control panel to
drive the wire pulls.  The six girdle sliders are written by the rhythm loop
every step, so dragging them has no effect.

## Motion acceptance test

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\test_cable_motions.py
```

Reports are written to `mujoco_arm/exports/cable_drive/`.

## Where guides may sit

## Control is coordinated, not one cable per joint

An arm is never driven by one muscle at a time, so the model must not be tuned
so that each cable acts on a single joint.  What matters is that the set of
sliders can produce the motions, which is what `test_cable_motions.py` checks
(it pulls the channels of one muscle group together) and what the viewer
exposes: nine independent `pull_<cable>` sliders plus six temporary girdle
sliders that the rhythm loop drives.

`optimize_guide_points.py` scores a guide on the worst signed moment arm about
its own joint over the working range.  `--coupling-weight` (default 0.5) then
only breaks ties: it keeps the candidates whose target moment arm is within
that fraction of the best available and picks the least coupled of those.  It
does not try to make a cable single-joint, and setting it to 0 scores the
target moment arm alone.

Passive joint stiffness stands in for the capsule and ligaments until real
ligaments are added (`defaults.joint_stiffness_nm_per_rad` in the manifest).
Without it the humeral twist is completely free, so any wire's twist moment
runs the shoulder to its +/-90 deg limit within half a second and the other
axes are carried around with it.

`test_slider_authority.py` sweeps every slider on its own and reports the joint
each one moves most, which is the quickest way to see whether a channel is
connected at all.

`optimize_guide_points.py` only samples the *surface* of an allowed bone: every
candidate is a triangle centroid pushed out along its normal by
`--clearance-mm` (default 6 mm), so the guide is a fixture that could actually
be bolted to a bone rather than a point floating off the arm.  The earlier
version searched a free box around the joint and produced fixtures 45-130 mm
off the bone surface, which no bracket can reach.

`--shoulder-origin` picks which bones the shoulder muscles may use:

- `torso` (default) bolts the fixture to the trunk frame.  The wire then spans
  only from the trunk to the humerus.
- `girdle` uses the anatomical origin bones (clavicle, scapula).
- `both` allows either and takes the better score.

With the girdle origins the wire loads the shoulder girdle far harder than the
arm: at 400 N the deltoid anterior develops 56 N.m about the clavicle axial
rotation but only 13 N.m about shoulder flexion, and the temporary rhythm
actuators are the only thing holding the girdle.  Use `girdle` again once the
girdle has real muscles.

`diagnose_pull.py` prints every skeleton angle while one cable is pulled, which
is how the shoulder failure was found: the wire's 3.7 N.m about the humeral
twist axis is unopposed, so the humerus twists to its +/-90 deg limit, the
flexion and abduction axes are carried around by that twist, and the shoulder
flexion stops at a few degrees even though the cable is at full force.

## Bone scaling

`exports/bone_scale.json` holds a uniform scale per bone.  The export applies
it to the part mesh and to every part-level coordinate system, about the part
origin.  The hand-drawn humerus and ulna place their part origin at the
proximal joint centre (`GH_J` and `ELBOW_J`), so scaling about that origin
keeps the joint centre fixed and lengthens the bone distally: the humerus went
from 220 mm to 264 mm at `x1.2` while the shoulder and elbow frames stayed put.

Set a bone back to 1.0 once the identical scale feature exists inside the
SolidWorks part, so the CAD and the simulation stay in step.  Re-run the
one-click batch afterwards; `check_bone_layout.py` re-verifies that every frame
still sits on its bone and that the joint gaps stay small.

## Guide optimization

After changing the five bone STLs or the joint frames, re-optimize the fixed
wire guides with:

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\optimize_guide_points.py
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\mujoco_arm\build_cable_drive_model.py
```

The optimizer writes explicit `guides` entries back into
`config/cable_drive_manifest.json`. Once those guides are frozen, the
one-click batch file uses them directly.

Current real-model status: the model compiles, all eight cable channels have
nonzero moment arms, and the one-way tension check passes. The scripted motion
test still reports partial shoulder range; the final guide-fixture tuning is
the remaining optimization work.

The builder deliberately fails if `joint_frames.json`, the five global bone
STLs or the cable manifest are missing. There is no fallback to the old CAD
model.
