# BodyParts3D arm: bones, real muscles, cables and printed joints

This is the second cable-drive pipeline.  It does not use the hand-drawn
SolidWorks skeleton at all: it takes the right arm straight out of
`D:\XJTLU\isa_BP3D_4.0_obj_99`, which is a body-wide model where every bone
and every muscle already sit in one coordinate frame and already articulate.

## One-click launch

- `open_bp3d_arm.bat` opens the model that is already built, in about a second.
- `run_bp3d_arm.bat` rebuilds everything, checks it, then opens the viewer.  It
  takes about a minute because the muscle centre lines are sliced out of the
  BodyParts3D meshes every time.

In the viewer every `pull_<wire>` slider is one wire, in metres.  Positive is a
pull; negative releases the wire.  The tendons are one-sided, so a released
wire carries nothing and cannot act as a spring against the others.

## Why the BodyParts3D skeleton

The SolidWorks route needed a coordinate system per joint, exported part
transforms, a bone scale file and a layout checker, and it still produced a
skeleton whose joints had to be trusted.  The BodyParts3D bones are already
assembled.  Measured gaps at the joints are 0.60 mm (sternoclavicular),
1.39 mm (acromioclavicular), 0.64 mm (glenohumeral), 1.06 mm (elbow) and
2.32 mm (proximal radioulnar), and the humerus is 307 mm long, so the model is
life size with no scaling step.

## Pipeline

```text
collect_bp3d_arm.py             bones, joint frames, collision proxies
build_bp3d_muscle_manifest.py   fused muscles and wire centre lines
build_bp3d_ligament_manifest.py ligament strips and joint ranges
build_bp3d_arm_model.py         arm_bp3d.xml
check_bp3d_arm.py               structure, moment arms, neutral pose
test_bp3d_motions.py            scripted cable motions
run_bp3d_arm.py                 viewer
```

### 1. Bones and joints

`collect_bp3d_arm.py` resolves element ids from `isa_element_parts.txt` by
concept name, so nothing is hard coded.  The thorax is fused from the right
1-12 ribs, the sternum, the cervical vertebrae and the thoracic vertebrae
T1-T12, 31 elements in total.

Joint centres are measured, not authored:

- glenohumeral: least-squares sphere through the humeral head, rms 2.92 mm,
  head radius 22.1 mm
- elbow: least-squares cylinder through the distal humerus, rms 3.45 mm,
  radius 11.8 mm
- sternoclavicular: midpoint of the closest clavicle / sternum vertex pair
- acromioclavicular: midpoint of the closest clavicle / scapula vertex pair
- proximal radioulnar: least-squares sphere through the radial head

The elbow axis is fitted at three band widths and averaged; the spread between
them is 6.7 degrees and that is the quality number written into the contract.
A single fit at 40 mm picks up the humeral shaft and returns a vertical axis,
which is why the multi-band average exists.

Bone mass comes from anatomical dry-bone values rather than from mesh volume
times cortical density: these meshes describe the outer surface of a whole
bone, and that over-estimates a real porous bone two to three times.  Inertia
is computed exactly from the mesh with the polyhedron second-moment sum and
written into the XML, so it does not depend on MuJoCo's own mesh inertia.

Anatomically, the joint axes are defined so that a positive rotation is the
named action.  For the right arm that required flipping three of them: turning
about the anterior direction raises the right arm medially, so abduction,
clavicle elevation and scapular upward rotation use the opposite sign, and
internal rotation uses the humeral axis pointing proximally.

### 2. Muscles and wires

Twenty-six BodyParts3D muscle elements are fused into 14 drive channels.  The
element list per channel is `bp3d_common.CHANNELS`.

The wire runs along the centre of the muscle volume, because the printed TPU
X-lattice wraps the wire.  The centre line comes from slicing the muscle
perpendicular to its principal axis and taking the centroid of each
cross-section.  A channel that needs more wires than it has muscles gets extra
parallel fibres offset sideways inside the largest belly, which is how several
Bowden cables sit side by side in one printed muscle.

Two corrections are applied automatically and both are reported in the
manifest:

- Role assignment.  An end is the origin when it is nearer the channel's
  declared origin bone.  A pure nearest-bone rule bolts the subscapularis to
  the ribs, because it lies only millimetres from them.
- Direction filtering.  A wire is kept only when `dL/dq` on the channel's
  target axis has the sign that makes a pull produce the named action.  The
  fibres of a wide fan such as the serratus anterior can end up on both sides
  of the axis and cancel each other; one wire was dropped this way.

### 3. Joint pulleys

A straight wire from a muscle belly to its insertion passes through the joint
centre once the joint bends far enough, and its moment arm changes sign: past
that angle the muscle pushes the joint the wrong way.  In the first build the
elbow wires crossed zero at 37 and 75 degrees, and the elbow refused to bend
past 50 degrees no matter how hard it was pulled.

Every wire therefore gets a wrapping point bolted to the joint's parent bone
at the joint cylinder, placed by an angular sweep that maximises the worst
moment arm over the joint's whole range with the correct sign.  After that
change every wire has a sign-stable moment arm: `dL/dq` does not cross zero
anywhere between the joint's limits, and the elbow reaches 147 degrees.

### 4. Ligaments: a 2D pattern projected onto the joint cylinder

The source framework builds a joint by unfolding the joint cylinder into a
plane, laying a 2D branch network on it and projecting it back.  A branch
behaves as a folded member: it unfolds freely and only becomes a stiff spring
once straight.  That is reproduced one dimension lower.

- Candidate branches are placed at six angles around each anatomical axis with
  the two ends offset by 75 degrees, so the branch is oblique to the axis the
  way a real collateral ligament is.  Without that offset the branch lies
  parallel to the axis and carries no moment at all.
- Each end is a real bone surface vertex inside a band around the joint
  cylinder, chosen by matching the angle and the design radius.
- A wrap point on the cylinder keeps the branch outside the bones.
- `springlength` is set so the branch goes taut at 86% or 100% of the joint's
  design limit, so the last degrees of every motion stiffen gradually instead
  of hitting a hard stop.
- Only branches that actually lengthen in a direction are kept.

81 strips are generated.  The hard joint range is set 5% outside the ligament
limit, so widening it further changes nothing: the ligaments are the real
constraint, which is the property the plan asks for.

Strip stiffness uses the printed TPU95A modulus of 26 MPa rather than the
folded-pattern modulus, because the folding is already represented by the free
length.  This is the one number to replace with a measurement from a printed
sample.

### 5. Cables and sliders

Each wire is a spatial tendon anchored on the insertion bone, wrapped over the
origin bone and finished on a slider.  The slider rides on the tail bone, not
on the world, which is the Bowden-tube case: the tube end is clamped to the
limb, so the wire's outlet moves with the body.

This matters more than it sounds.  With a world-fixed slider every wire gets a
lever arm of 100-250 mm/rad about the sternoclavicular joint, and pulling any
channel simply pivots the whole arm from the base.  The first world-fixed build
gave the elbow flexors a 138 mm/rad moment about the clavicle against 16 mm/rad
about the elbow, and the elbow stopped at 50 degrees while the girdle swung
instead.

MuJoCo's tendon spring is two-sided, so a slack wire would push.  Every wire
and every ligament carries an `assist_` general actuator that cancels the
compressive half.  The sign convention was verified on a minimal model: with a
tendon shorter than its `springlength` the spring produces a *positive*
generalized force, so the correction is `min(0, k (L - rest))` and the actuator
range has to include negative values.  The legacy helper in
`cable_drive_common.py` uses the opposite sign and was left alone so the older
model keeps behaving the way it was signed off.

## Measured behaviour

`check_bp3d_arm.py` on the current model:

```text
bones 6, skeleton joints 5, wires 22, ligaments 81, sliders 22
joint gaps 0.60 ... 2.32 mm, collars inside the bone band
neutral pose finite, no taut ligament
rest pose  sternoclavicular 7.8  acromioclavicular 0.8
           glenohumeral 13.7  elbow 10.3 deg
every wire has a sign-stable moment arm across its range
```

`test_bp3d_motions.py`, gravity on, four second ramp, antagonist released:

```text
elbow flexion              147.3   target 120    pass
shoulder abduction          71.4   target  60    pass
shoulder flexion            63.9   target  90    below
shoulder external rotation -24.4   target -30    pass
shoulder internal rotation  61.0   target  30    pass
scapular upward rotation    52.7   target  25    pass
clavicle elevation          -3.6   target  15    below
elbow extension return      -0.6   target <=20  pass
```

The motion test applies a co-contraction fraction to the antagonist rotator
cuff during the shoulder motions.  That is not a trick to make the numbers look
better: a ball joint driven by one cable takes the shortest path and turns the
motion into a twist, and holding the humeral head in the glenoid plane is what
the cuff does.  Leaving it out is what a single-cable test measures, and it is
not how the device is commanded.

## Control interface

```text
pull_<wire>   slider position in metres, positive pulls, negative releases
```

Wires are grouped by channel.  A channel with two wires is one muscle with two
Bowden cables and both should be moved together.

## Known limitations

- Shoulder flexion reaches 64 of 90 degrees.  The deltoid anterior and
  coracobrachialis moment arms are 15-16 mm/rad, and the ball joint still
  spends part of the pull on its other axes.  More travel, or a stronger
  anterior deltoid group, is the next thing to try.
- Clavicle elevation goes slightly negative.  Pulling the trapezius rotates
  the clavicle about its other axes as well, so the elevation component alone
  stays small.  This needs the two clavicle channels commanded together.
- Scapulothoracic contact is off.  The convex hull of the scapula fills in the
  concavity that lets the real bone sit on the rib cage, so the hull intrudes
  up to 14 mm into the rib proxies at the neutral pose.  Contact is enabled
  only for the pairs whose proxies are faithful: the humeral head sphere
  against the scapula hull, and the ulna against the radius.  A concave
  decomposition of the scapula and the clavicle is what unblocks the sliding
  interface.
- The muscle bodies are not simulated yet.  The printed X-lattice is a visual
  asset in this version: the wire runs down its centre but the lattice does
  not carry load.  The next step is to replace the visual mesh with a `flex`
  body and couple the tendon strain to the lattice nodes.
- The ligament free lengths are geometry, not measurement.  Everything follows
  from the joint cylinder, the strip angles and the design limits.  Print one
  strip and the free-length model can be replaced with a real curve.
- The joint collar stiffness is a chosen value, in `DEFAULT_JOINT_STIFFNESS`
  in `build_bp3d_arm_model.py`.  It stands in for the printed TPU collar, which
  is what holds the limb between cable commands.  The girdle is stiff at
  8 N.m/rad so the arm droops only a few degrees with every wire slack, while
  the humerus and elbow are left free at 2.5 and 1.5 N.m/rad so the collar
  does not eat the cable torque.

## Parameters carried from the source framework

```text
lattice, deformable belly   (x, y, z, p) = (14, 12, 30, 2) mm
lattice, stiff ends         (x, y, z, p) = ( 6,  8, 25, 2) mm
muscle cell compressive stiffness    1.3 N/mm
tendon cell tensile stiffness       15.0 N/mm
printed strip                8 x 2 mm, printed TPU95A, 26 MPa
Bowden outlet                bolted to the moving bone
```

These are the values the source framework measured for printed TPU95A cells.
They are recorded per channel in `config/bp3d_muscle_manifest.json` so the
Grasshopper generator and the simulation stay on the same numbers.

## Files

```text
bp3d/*.obj                        bones and collision hulls, metres, global frame
bp3d/muscles/<channel>_<n>.obj    muscle volume behind each wire
config/bp3d_arm.json              bones, joint frames, collision proxies
config/bp3d_muscle_manifest.json  14 channels, 22 wires, centre lines
config/bp3d_ligament_manifest.json 81 strips
config/bp3d_joint_rom.json        per joint, axis and direction strip report
arm_bp3d.xml                      the MuJoCo model
exports/bp3d/*.json               check, motion test and wire reports
```

The older SolidWorks pipeline (`arm_cable_drive.xml`, `open_arm.bat`,
`run_cable_drive.bat`) is untouched and still runs.
