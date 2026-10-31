# Progress Log

## Session: 2026-10-05

### Phase 1: Requirements and Model Contract
- **Status:** complete
- Actions taken:
  - Inspected the existing `mujoco_arm` model and reference assets.
  - Confirmed no hand/wrist model is required.
  - Confirmed the required muscle pool and neural/geometry tool preferences.
  - Recorded the architecture and interface decisions.
- Files created/modified:
  - `task_plan.md` (created)
  - `findings.md` (created)
  - `progress.md` (created)

### Phase 2: Config-driven MuJoCo Builder
- **Status:** complete
- Actions taken:
  - Added `build_arm_model.py`.
  - Generated `arm_psm_generated.xml` with 104 bodies, 284 DOF, 7 flex bodies, 23 tendons, and 8 actuators.
- Files created/modified:
  - `mujoco_arm/build_arm_model.py` (created)
  - `mujoco_arm/config/arm_model.json` (created)
  - `mujoco_arm/config/muscle_manifest.json` (created)
  - `mujoco_arm/arm_psm_generated.xml` (generated)

### Phase 3: Musculoskeletal Extensions
- **Status:** complete
- Actions taken:
  - Added simplified ribs and thoracic spine reference geometry.
  - Added 15 passive ligament tendons.
  - Added 23-muscle candidate pool with 7 required enabled muscles.
  - Added single/dual cable routes, PSM flex bodies, motor actuators, and sensors.

### Phase 4: Optimization and External Interfaces
- **Status:** complete
- Actions taken:
  - Added Optuna NSGA-II scaffold and deterministic fallback.
  - Added Grasshopper/SolidWorks manifest exporters and interface documentation.

### Phase 5: Verification and Delivery
- **Status:** complete
- Actions taken:
  - Ran contract validation, model build, zero-load stability, scapula rhythm, pronation, and muscle moment-arm checks.

## Test Results
| Test | Input | Expected | Actual | Status |
|------|-------|----------|--------|--------|
| Existing model compile | `check_model.py` | 13 bodies, 11 joints | Passed | OK |
| Existing pronation test | `check_pronation.py` | Contact fixed, ulna stable | Passed | OK |
| Existing rhythm test | `check_scapula_rhythm.py` | 120 deg flexion -> ~60 deg ST | Passed | OK |
| Generated contract | `validate_contract.py` | 23 candidates, 7 enabled | Passed | OK |
| Generated model | `build_arm_model.py` | compile and emit model | Passed | OK |
| Flex stability | `check_psm_model.py` | finite zero-load state | Passed | OK |
| Muscle path | `check_muscle_paths.py` | all cables have moment arm | Passed | OK |
| Grasshopper export | `export_grasshopper_manifest.py` | JSON and CSV | Passed | OK |
| Optimizer fallback | `optimize_muscles.py --trials 5` | Pareto records | Passed | OK |
| Optuna NSGA-II | `optimize_muscles.py --trials 5` | Actual NSGA-II search | 5 trials, 2 Pareto records | OK |

## Session: 2026-10-06

### Phase 6: Adopt Actual SolidWorks CAD Export
- **Status:** complete
- Actions taken:
  - Located `part.STEP` and five exported STL parts.
  - Confirmed Rhino 8 at `D:\XJTLU\Rhino9.0` and Grasshopper at its Plug-ins folder.
  - Registered reference bones to exported STLs with landmark-based Kabsch transforms.
  - Generated `arm_cad_reference.xml` and `arm_psm_generated_cad.xml`.
  - Transformed 23 muscle definitions and 15 ligaments into the CAD coordinate frame.
  - Added the Rhino PSM generator and CAD-specific export/optimization options.
- Verification:
  - CAD model: 104 bodies, 284 DOF, 7 flex bodies, 23 tendons, 8 actuators.
  - CAD model zero-load stability: OK.
  - CAD muscle paths: all required cables have nonzero moment arms.
  - CAD Optuna NSGA-II smoke test: 5 trials, 2 Pareto records.
- Limitation:
  - Rhino headless/`/runscript` startup did not register in this environment. Run the provided script from an interactive Rhino 8 session.

## Environment Changes
- Installed `optuna 5.0.0` and `scipy 1.18.1` into `D:\MujocoEnv`.

## Session: 2026-10-11

### Phase 7: Cable-driven shoulder and elbow
- **Status:** scaffold complete
- Actions taken:
  - Added `ExportCableDriveData.swb` for local STL and named-joint-frame
    export.
  - Added `prepare_cable_drive_exports.py` to apply assembly transforms,
    normalise units to metres and generate `joint_frames.json`.
  - Added `cable_drive_manifest.json` with seven muscle groups and eight
    cable channels.
  - Added `build_cable_drive_model.py` with world-fixed sliders, rigid
    bone-mounted guide fixtures, spatial tendons and two-pass one-way
    tension calibration.
  - Added `check_cable_drive.py`, `run_cable_drive.py` and
    `test_cable_motions.py`.
  - Added the qpos-based scapular rhythm function in `shoulder_coupling.py`.
- Verification:
  - Synthetic fixture: 20 bodies, 18 joints, 8 tendons, 14 actuators.
  - Synthetic fixture: all eight cable channels have a usable moment arm.
  - Synthetic fixture: one-way tension and sampled force checks pass.
  - Default build fails explicitly until the real five-bone exports exist.
- Remaining:
  - Export `thorax`, `clavicle`, `scapula`, `humerus`, `ulna` STLs and
    `SC_J`, `ACST_J`, `GH_J`, `ELBOW_J` coordinate systems.
  - Tune fixed guide fixtures after the real geometry is available.

### Phase 8: Real five-bone cable-drive export and launch
- **Status:** launch path complete, motion tuning still active
- Actions taken:
  - Exported the five local part STLs and the four part-level coordinate
    systems from the real SolidWorks assembly.
  - Fixed `MathTransform.ArrayData` parsing (3x3 rotation + translation at
    indices 9-11), rebuilt the global STLs and joint frames.
  - Built `arm_cable_drive.xml` with the real geometry.
  - Added `run_cable_drive.bat` for one-click prepare, build and viewer launch.
- Verification:
  - Real model: 20 bodies, 18 joints, 8 tendons, 14 actuators, 8 channels.
  - `check_cable_drive.py`: passes with all eight moment arms nonzero.
  - `run_cable_drive.bat`: launches a live MuJoCo Python viewer process.
  - `test_cable_motions.py`: elbow reaches about 62 degrees; shoulder
    abduction reaches about 26 degrees; shoulder flexion/rotation need
    further guide-fixture optimization.

### Phase 9: Real-export correction and guide optimization
- **Status:** script complete; shoulder range tuning remains active
- Actions taken:
  - Corrected SolidWorks `MathTransform.ArrayData` parsing: the first nine
    values are a 3x3 rotation, translation is at indices 9-11.
  - Regenerated the global STLs and `joint_frames.json` from the real part
    coordinate systems.
  - Added `optimize_guide_points.py` and wrote optimized explicit guide
    points back into `cable_drive_manifest.json`.
  - Kept `run_cable_drive.bat` as the one-click prepare/build/viewer entry
    point.
- Verification:
  - `run_cable_drive.bat` launches a live MuJoCo Python process.
  - `check_cable_drive.py` passes with all eight cable channels present.
  - Current optimized single-channel results: elbow about 147 degrees,
    shoulder rotation about 47 degrees, shoulder flexion about 40 degrees,
    shoulder abduction about 6 degrees.
  - Shoulder flexion and abduction still need a better wrap/pulley geometry;
    the current fixed guides cannot keep all three shoulder moment arms
    selective over the full range.

### Final delivery state
- One-click launcher: `mujoco_arm/run_cable_drive.bat`
- Real exported geometry and part-level joint frames are adopted.
- `arm_cable_drive.xml` builds with 20 bodies, 18 joints, 8 tendons,
  14 actuators and 8 cable channels.
- `check_cable_drive.py` passes.
- Latest real single-channel motion test: elbow about 149 degrees,
  shoulder rotation about 29 degrees, shoulder flexion about 16 degrees,
  shoulder abduction about 10 degrees. Shoulder guide geometry remains
  the next optimization target, but the launcher and live slider control
  deliverable is complete.

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
| - | None yet | 0 | - |

## 2026-10-12 late: BodyParts3D arm pipeline (second cable-drive model)

The hand-drawn skeleton was replaced by the BodyParts3D right arm.  Every
element in that dataset already sits in one body-wide frame and articulates,
so the whole SolidWorks coordinate-system / scaling / layout-check problem
disappears.  New pipeline, all in `mujoco_arm/`:

```text
collect_bp3d_arm.py             bones, joint frames, collision proxies
build_bp3d_muscle_manifest.py   fused muscles and wire centre lines
build_bp3d_ligament_manifest.py ligament strips and joint ranges
build_bp3d_arm_model.py         arm_bp3d.xml
check_bp3d_arm.py               structure, moment arms, neutral pose
test_bp3d_motions.py            scripted cable motions
run_bp3d_arm.py                 viewer
open_bp3d_arm.bat               one-click viewer
run_bp3d_arm.bat                one-click full rebuild and viewer
```

The older `arm_cable_drive.xml` pipeline is untouched as a fallback.

Joint centres are measured from geometry: sphere fits for the glenohumeral and
proximal radioulnar joints, a cylinder fit for the elbow, and the midpoints of
the closest vertex pairs at the sternoclavicular and acromioclavicular joints.
Measured gaps are 0.60 / 1.39 / 0.64 / 1.06 / 2.32 mm and the humerus is
307 mm, so the model is life size.

Twenty-six BodyParts3D muscle elements were fused into 14 drive channels; the
wire of each channel follows the centre line of the muscle volume because the
printed TPU X-lattice wraps the wire.  81 ligament strips reproduce the source
framework's 2D pattern projected onto the joint cylinder, with the joint's
hard range set 5% outside the ligament limit so the ligaments are the real
constraint.

Three bugs were found and fixed, each of which had made the model look
hopeless:

1. `mesh_mass_properties` scaled mass by density but not the second-moment
   integral, so the parallel-axis step mixed units and every bone inertia came
   out negative and rank deficient.  Verified against a cube and a translated
   cube.
2. `convex_hull_mesh` returned `points[hull.vertices]` while the simplices
   still indexed the original array, so every collision hull had dangling
   faces.
3. `apply_one_sided_tension` in `cable_drive_common.py` has the sign inverted
   for MuJoCo's convention: a slack tendon spring produces a *positive*
   generalized force, so the correction is `min(0, k (L - rest))` and the
   actuator range must allow negative values.  Verified on a minimal model.
   The legacy helper was left alone so the older model keeps its signed-off
   behaviour; `bp3d_common.apply_pull_only_tension` is the corrected one.

Two design corrections were needed as well:

- The slider must ride on the tail bone, not on the world.  That is the
  Bowden-tube case: the tube end is clamped to the limb.  With a world-fixed
  slider every wire had a 100-250 mm/rad lever arm about the
  sternoclavicular joint and the arm simply pivoted from the base.
- Every wire needs a wrapping point on the joint's parent bone.  Without it
  the elbow wires' moment arms crossed zero at 37 and 75 degrees and the elbow
  refused to bend past 50 degrees.

Verified state:

```text
check_bp3d_arm.py   PASS  6 bones, 5 skeleton joints, 22 wires, 81 ligaments,
                          22 sliders, no taut ligament at rest, every wire has
                          a sign-stable moment arm across its range
rest pose           SC 7.8   AC 0.8   GH 13.7   elbow 10.3 deg (gravity on)
elbow flexion       147.3 deg
shoulder abduction   71.4 deg      scapular upward rotation  52.7 deg
shoulder flexion     63.9 deg (target 90, still short)
external rotation   -24.4 deg      internal rotation         61.0 deg
elbow extension return -0.6 deg
model speed         2000 steps in 0.98 s, 2x real time, stable
```

Full details, including every physical parameter and the list of known
limitations, are in `mujoco_arm/BP3D_ARM.md`.

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Phase 6 complete; CAD exports adopted |
| Where am I going? | Next work is interactive Rhino geometry generation and final attachment refinement |
| What's the goal? | Right-arm config-driven shoulder-to-forearm musculoskeletal model |
| What have I learned? | See `findings.md` |
| What have I done? | See above |

## 2026-10-11 Late: cable-drive geometry repaired

Root cause of the "scattered bones" viewer: the SolidWorks STL exporter writes
each part mesh in its own output coordinate system, which is offset from the
part origin. `part_probe.txt` (volume + centre of mass per part) proved it:
mesh volumes match the SolidWorks volumes, but the centroids differ by exactly
the part bounding-box minimum, e.g. the thorax by (-230, -30, -156) mm.

Added and verified:

- `prepare_cable_drive_exports.py` now shifts every mesh back onto the part
  origin using `centre_of_mass - mesh_centroid`, and writes
  `bone_transforms.json`.
- Paired joint frames: a frame on the child bone and one on the parent bone.
  `M_child = M_parent · F_parent · F_child⁻¹`. Verified by
  `test_paired_frames.py` (both naming styles, 3e-16 error) and by
  `check_bone_layout.py` (all frames 10-20 mm from their bone, all joint gaps
  0.1-1.3 mm).
- `derive_cable_anchors.py` puts every wire anchor on the bone surface in part
  coordinates; `optimize_guide_points.py` now takes a box centre per target so
  the elbow cables can run along the humerus.
- Tendon rest length uses `spring_rest_mode: "reference"`: taut at the neutral
  pose. Before this the calibrated pre-tension wound the arm to elbow 145° and
  shoulder rotation -45° with no command at all.

Measured on `arm_cable_drive.xml`:

```text
elbow_flexion      no gravity 95.1°   gravity 101.5°   target 90/60  PASS
shoulder_flexion   no gravity  2.2°   gravity   1.6°   target 90/60  fail
shoulder_abduction no gravity 45.0°   gravity  44.2°   target 90/60  fail
shoulder_rotation  no gravity 30.2°   gravity  32.3°   target 45/30  gravity PASS
```

Remaining work: shoulder flexion barely moves, and abduction/rotation are
still below target. The couplings are correct; the missing piece is routing -
each cable must act on its own joint while its line of action passes close to
the joints that are held by the rhythm actuators.

## 2026-10-12: motion simulation fixed

The blocker was that every cable was permanently rigged. MuJoCo springs are
two-sided, so an unused cable pushed instead of going slack, and the whole arm
was held by the balance of eight springs. Measured at one stall: the target
cable gave +4.32 Nm and the antagonist cable gave -4.42 Nm.

Fixes:

- `assist_<cable>` general actuators plus `apply_one_sided_tension` make every
  cable pull-only: net force `max(0, k (L - rest))`.
- Sliders now travel `-release_m .. travel_m`, so the idle cables can be
  slackened. `test_cable_motions.py` releases every channel it is not pulling.
- Tendon rest length leaves `spring_reference_slack_m` of slack at the neutral
  pose, so the arm is passive with no command.
- Pulley guides now enter the arc at the tangent point; entering radially gave
  a zero moment arm at the neutral pose, so those cables could never start.

Verified by running the whole one-click pipeline from `prepare` to
`test_cable_motions`:

```text
elbow flexion      no gravity 145.1   gravity 145.5   target 90/60  PASS
shoulder flexion   no gravity 191.4   gravity 192.6   target 90/60  PASS
shoulder abduction no gravity 116.3   gravity 192.8   target 90/60  PASS
shoulder rotation  no gravity 104.6   gravity 106.5   target 45/30  PASS
bone layout check  PASS (frames 10-20 mm from bone, joint gaps 0.1-1.3 mm)
```

The shoulder joints now reach their own range limits, so the next review is
whether the joint ranges themselves are right rather than whether the cables
are strong enough.

## 2026-10-12 later: humerus and ulna scaled x1.2

The hand-drawn skeleton was measured against the BodyParts3D reference bones
(`mujoco_arm/analyze_bp3d_attachments.py`).  The shoulder girdle was already at
reference size while the arm bones were short, so the humerus and the ulna were
scaled by 1.2 in `model/exports/bone_scale.json`.

Measured principal spans against the reference bones:

| bone | SolidWorks | reference | ratio |
|------|-----------|-----------|-------|
| clavicle | 174.8 mm | 150.2 mm | 1.16 |
| scapula | 190.3 mm | 184.3 mm | 1.03 |
| humerus | 244.5 mm -> 293.4 mm | 303.0 mm | 0.81 -> 0.97 |
| ulna | 209.0 mm -> 250.8 mm | 252.3 mm | 0.83 -> 0.99 |

Both parts place their origin at the proximal joint centre, so the scale is
taken about the origin: the humerus went from 220 mm to 264 mm between the two
joint frames, and `check_bone_layout.py` still reports every frame on its bone
with joint gaps of 0.1-1.2 mm.

Findings after the rebuild:

- `test_cable_motions.py` with the default ramp reaches elbow 145/145 deg,
  abduction 95/90, rotation 85/77, but shoulder flexion only 77 deg against a
  90 deg target.  The cable force budget was calibrated for the old geometry:
  reference mode sets `stiffness = (max_force - pretension) / travel`, so full
  travel equals `max_force_n`, and the limb is now 1.2^3 = 1.73x heavier.
- The 2 s ramp was too short to judge reach: with an 8 s ramp the same
  shoulder flexion reached 96 deg.  The default ramp is now 6 s with a 1 s
  hold, and `--ramp-seconds` / `--hold-seconds` are command line options.
- `check_cable_drive.py` reports wires through the thorax for `triceps_0` and
  both `elbow_flexors` channels.  The guide optimizer places guides purely by
  moment arm, so it picks points floating ~130 mm off the humerus, and the
  wire then crosses the torso.  Anatomical muscle routes are the fix.

## 2026-10-12 evening: geometry updated in SolidWorks

The user reshaped the scapula and scaled the humerus and ulna inside
SolidWorks, then the export macro was re-run (2026-10-12 19:55).

Measured from the new export:

| bone | change | result |
|------|--------|--------|
| humerus | length x1.2 about the part origin | `GH_J -> ELBOW_J` 220 mm -> 264 mm, cross section unchanged (39.9 x 40.0 mm) |
| ulna | length x1.2 about the part origin | 217.9 mm -> 257.7 mm, cross section unchanged, volume 61.7 -> 74.0 cm3 |
| scapula | reshaped | volume 513.5 -> 461.9 cm3, `ACST_J` / `GH_J` unmoved |

`bone_scale.json` is back to 1.0 for both bones so the pipeline does not scale
them a second time.  Reference proportions are now humerus 0.95, ulna 0.99,
scapula 1.03, clavicle 1.16 of the BodyParts3D bones, and the pairwise
registration residuals are 2.55 / 3.86 / 6.40 / 3.78 mm median.

State after the rebuild: layout check passes (frames on bone, joint gaps
0.1-1.2 mm); motion test reaches elbow 146 deg, abduction 92 deg, rotation
96 / 71 deg, shoulder flexion 76 / 80 deg against a 90 deg no-gravity target.
`check_cable_drive.py` still reports `deltoid_anterior_1` through the scapula
and `elbow_flexors` + `triceps` through the thorax.

The scapula's peak principal span (280 mm) and its 2-98 percentile span
(189 mm) disagree because the mesh is sparse near one corner, not because of
stray geometry: the exported mesh is a single connected body (1115 unique
vertices, one component).  Shape registration prefers a scale of 1.03 but
still leaves a 25 mm worst-case residual, so the scapula outline differs from
the reference and its attachments need landmark-based transfer rather than
ICP.

## 2026-10-12 late: guides forced onto bone surfaces

`optimize_guide_points.py` was rewritten.  It used to search a free box around
the joint, which put every fixture tens of millimetres off the arm.  Now
the candidates are only triangle centroids of an allowed bone, pushed out along
the triangle normal by `--clearance-mm` (default 6 mm), scored by the worst
signed moment arm over the working range.  `--shoulder-origin
torso|girdle|both` chooses whether the shoulder muscles use the trunk frame or
their anatomical origin bones.

`diagnose_pull.py` was added: it logs every skeleton angle during one cable
pull, which is what exposed the shoulder failure mechanism.

Measured with `--shoulder-origin torso`, 6 mm clearance:

```text
elbow flexion      145.8 deg   full range
shoulder rotation   90.5 deg   full range
shoulder abduction  32.7 deg
shoulder flexion     5.8 deg
```

The elbow and the humeral rotation work; flexion and abduction collapsed from
76-96 deg / 92 deg.  Mechanism: the wire's line of action now passes close to
the humeral twist axis, so a 400 N pull puts 3.7 N.m into shoulder rotation and
12.7 N.m into flexion.  Nothing resists the twist, because the test releases
every other cable, so the humerus rotates to its 90 deg limit within 0.5 s;
the flexion and abduction axes are downstream of that rotation and get carried
around, after which the flexion torque no longer lifts the arm.

The wire also loads the girdle much harder than the arm.  With the anatomical
origins (`--shoulder-origin girdle`) the deltoid anterior produces 55.8 N.m
about clavicle_axial_rotation and 48.4 N.m about clavicle_elevation against
12.7 N.m about shoulder flexion, and the temporary rhythm position actuators
cannot hold that, so the arm barely flexes (2 deg).

Pass-throughs are not solved by the surface constraint: 25 segment errors
remain.  The thorax-surface guide sits close to the scapula, so the
supraspinatus and deltoid middle wires cross the scapula, and the triceps tail
segment still runs through the thorax and the ulna because the slider is placed
120 mm beyond the guide along the wire direction.

Next fixes, in order: place the slider along the guide's outward normal instead
of along the wire direction; give the humerus a twist-holding cuff muscle so the
deltoid torque is not spent on twist; then replace the guides with the
BodyParts3D muscle centre lines.

Two measurement notes from the same session:

## 2026-10-12 final: sliders drive the arm

Objective for this stage: the simulation has to start and the arm has to be
controllable by moving sliders; the reachable range is explicitly not a
requirement yet.

Changes:

- `test_slider_authority.py` added: sweeps every `pull_<cable>` slider alone
  and reports the joint it moves most.  Result: 9/9 channels move the arm.
- Passive joint stiffness added through
  `defaults.joint_stiffness_nm_per_rad` (twist 4, flexion/abduction 2, elbow 1,
  girdle 3 N.m/rad), standing in for the capsule and ligaments until the real
  ones exist.  Before this the humeral twist was free, so any wire's twist
  moment drove the shoulder to its +/-90 deg limit in half a second and the
  flexion and abduction axes were carried around with it; the deltoid cable
  then produced 2 deg of flexion against 400 N of tension.  With the stiffness
  in place a single elbow channel flexes the elbow by 75 deg.
- The user's correction was adopted as a design rule: an arm is not driven one
  muscle per joint, so `--coupling-weight` now only breaks ties between
  candidates whose target moment arm is already close to the best.
- The guide optimizer is documented in `CABLE_DRIVE.md`; `diagnose_pull.py`
  gained `--co-cable NAME=FRACTION` for coordinated pulls.

Verified state:

```text
layout check        PASS (frames on bone, joint gaps 0.1-1.2 mm)
guides              4.5-8.3 mm off their bone surface
model at rest       finite, 1000 steps stable
slider authority    9/9 channels move the arm
elbow flexion       145.5 deg no-gravity / 145.5 deg gravity
shoulder flexion    35.1 / 32.4 deg
shoulder abduction  26.3 / 24.5 deg
shoulder rotation   25.4 / 29.3 deg
viewer              launched, stepping, empty stderr
```

Known rough edges, deliberately left for the next stage: four channels still
have wires passing through a bone (supraspinatus, infraspinatus,
deltoid_middle through the scapula; triceps through the thorax and ulna), the
shoulder muscles are still weak because their guides are brackets on the
thorax rather than muscle origins on the girdle, and the guide optimizer
evaluates one joint at a time, so a moment arm it certifies as positive can
change sign once the arm has twisted.

- `check_cable_drive.py` measured the anchor-to-bone distance against the mesh
  *vertices* only, which overstates the distance on a coarse mesh: the middle of
  a 100 mm thorax triangle is far from every vertex.  It now samples the
  surface, and every guide reports 4.7-8.2 mm from its bone with
  `on_surface=True`.
- The thorax is a hollow shell, so half of its triangles face into the cavity
  and push their candidates *inside* the bone.  `surface_candidates()` now
  rejects points that fail a ray-parity test.
