# SolidWorks interface

The final CAD model must provide the following named coordinate systems:

- `SC_*`: sternoclavicular joint
- `AC_*`: acromioclavicular joint
- `GH_*`: glenohumeral joint
- `ELBOW_*`: elbow centre and axis
- `PRUJ_*`: proximal radioulnar joint
- `DRUJ_*`: distal radioulnar joint
- `WRIST_END_*`: forearm endpoint, with no hand or wrist DOF

For each enabled muscle, the CAD model should include:

- a conformal attachment pad on the origin and insertion bones
- a cable exit and central through-channel
- guide or pulley geometry at every listed via point
- a motor or linear-slide mounting interface at the root

The coordinates exported from SolidWorks must be expressed in metres in the
same rest pose as `arm_reference.xml`.
