# Findings

- Source STL and OBJ share shapes but differ by a common XY translation;
  calculate and verify this per bone rather than mixing coordinates.
- Source STL bounds and existing SolidWorks verification logs support mm units.
- Source has five anatomical bones and no thorax model.
- Earlier V2-S1 motion results and wire anchors do not transfer to A1.
- The user selected full arm architecture plus detailed elbow muscles first,
  native source dimensions, and no forearm pronation/supination.

- OBJ-to-assembly-STL translation approximately (282.33,148.417,0) mm for all
  five bones; individual differences are recorded, no independent recentering.
- Native source bone pairs have no solid penetration in neutral. GH minimum
  gap ~0.576 mm, humerus/ulna ~0.246 mm. Source bone form is retained unchanged.
- A collision-informed elbow axis adjustment of 2.457 mm / 2.206 degrees
  extends the sampled bone-only prefix from 41.75 to 145 degrees.
- Bone-only prefixes: GH flex120, abd118.75, int70/ext90; prescribed scapular
  up50; actual combined humerus direction168.08. No hardware acceptance.
- SC2/AC1 cannot reproduce the additional pure scapular translation/normal-axis
  rotation tasks within the tested pose tolerance beyond tiny prefixes.
- Two shaped lattice specimens are closed single material components and
  neutral bone-clear. This does not validate channel collapse/elongation.
- Float32 mesh extraction introduced micron-size slivers that invalidated
  OpenCascade sewing. Use Manifold.to_mesh64 for CAD input and recheck export.
