# Motion validation for `biomimetic_structure.SLDASM`

This check separates CAD mechanism validity from final musculoskeletal
performance.

## Current audit result (2026-10-10)

The read-only audit currently finds:

- 5 active components: `thorax`, `clavicle`, `scapula`, `humerus`, `ulna`
- no `radius`
- 1 concentric mate
- 0 reference axes
- 0 coordinate systems
- 0 features named for `SC_`, `AC_`, `ST_`, `GH_`, `ELBOW_`, `PRUJ_` or
  `DRUJ_`

This version is not ready for a biomimetic motion check. It can display the
bones in an assembly pose, but it cannot independently test scapulohumeral
rhythm, elbow motion and forearm pronation/supination.

The first required geometry additions are:

1. Add the `radius` part and position it in the neutral forearm pose.
2. Create named reference axes and coordinate systems for every joint.
3. Give the model one clear, independently testable motion chain before
   adding muscles.

Recommended reference names:

- `SC_ORIGIN`, `SC_AXIS_*`
- `AC_ORIGIN`, `AC_AXIS_*`
- `ST_ORIGIN`, `ST_AXIS_*`
- `GH_ORIGIN`, `GH_AXIS_*`
- `ELBOW_ORIGIN`, `ELBOW_AXIS`
- `PRUJ_ORIGIN`, `PRUJ_AXIS`
- `DRUJ_ORIGIN`, `DRUJ_AXIS`

## 1. Run the read-only SolidWorks audit

Open `biomimetic_structure.SLDASM` in SolidWorks and run:

`Tools > Macro > Run > D:\XJTLU\PSP\model\validation\InspectBiomimeticAssembly.swb`

The macro does not save the assembly. It writes:

- `assembly_audit.txt`
- `components.txt`
- `features.txt`
- `mates.txt`

## 2. Run the contract check

```powershell
D:\MujocoEnv\Scripts\python.exe D:\XJTLU\PSP\model\validation\validate_motion_contract.py
```

The result is written to `motion_contract_report.json`.

## 3. Check the mechanism in SolidWorks

Use `Evaluate > Interference Detection` at the neutral pose and at the joint
limits. Then use `Motion Study` to move only the intended joint and verify that
unintended bones remain fixed.

SolidWorks is the geometry and interference check. MuJoCo should be the
authoritative dynamic test because the shoulder chain and radioulnar coupling
are difficult to constrain correctly with SolidWorks mates alone. Avoid adding
so many mates that the model becomes over-defined; export the named axes,
centres and attachment frames and reproduce the intended degrees of freedom in
MuJoCo.

The intended kinematic chain is:

1. thorax -> clavicle: sternoclavicular motion
2. clavicle -> scapula: acromioclavicular motion
3. thorax/scapula interface: scapulothoracic sliding and rotation
4. scapula -> humerus: glenohumeral ball motion
5. humerus -> ulna: elbow hinge
6. ulna -> radius: proximal and distal radioulnar motion

The minimum motion acceptance tests are:

- shoulder flexion/abduction to approximately 180 degrees with a scapulohumeral
  rhythm close to 2:1
- elbow flexion/extension from 0 to approximately 145 degrees
- forearm pronation/supination close to +/-90 degrees, with the ulna stable
  and the radial head remaining at the capitular contact
- no component separation, unintended penetration or over-constraint
- key corresponding points in SolidWorks and MuJoCo within 1 mm and joint
  axes within 1 degree

## 4. Export for MuJoCo

The minimum export set is one active part for every bone shown above. If
`radius` is absent, forearm pronation and supination cannot be validated.

Export each part in its assembly pose, or export the individual parts plus
their component transforms, then register them into the same assembly frame.
The exact joint centres and axes should be exported as named reference points,
axes or coordinate systems using the prefixes `SC_`, `AC_`, `GH_`,
`ELBOW_`, `PRUJ_` and `DRUJ_`.
