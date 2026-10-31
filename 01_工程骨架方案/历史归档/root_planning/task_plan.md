# Task Plan: Right-arm shoulder-to-elbow X-lattice musculoskeletal model

## Goal
Implement a config-driven right-arm MuJoCo framework for the shoulder, upper arm, elbow,
and forearm with thoracic reference geometry, passive constraints, a muscle pool,
X-lattice flex bodies, cable routes, and an Optuna/Grasshopper integration scaffold.

## Current Phase
Phase 7 scaffold complete; waiting for the new five-bone exports

## Phases

### Phase 1: Requirements and Model Contract
- [x] Confirm right arm only
- [x] Exclude hand and wrist joints
- [x] Confirm SolidWorks plus MuJoCo workflow
- [x] Confirm MuJoCo flex, Rhino/Grasshopper, and Optuna choices
- [x] Define unified JSON schemas and placeholder geometry
- **Status:** complete

### Phase 2: Config-driven MuJoCo Builder
- [x] Implement `arm_model.json` and `muscle_manifest.json`
- [x] Implement Python builder for skeleton, joints, pectoral girdle, and forearm
- [x] Generate and compile a MuJoCo model from the JSON files
- **Status:** complete

### Phase 3: Musculoskeletal Extensions
- [x] Add thorax/rib/spine attachment geometry
- [x] Add passive tendon/equality constraints for ligaments
- [x] Add muscle pool, PSM shape classes, flex bodies, and cable routes
- [x] Add stable initial controls and sensors
- **Status:** complete

### Phase 4: Optimization and External Interfaces
- [x] Add Optuna NSGA-II scaffold with design variables and constraints
- [x] Add Grasshopper/SolidWorks export contract and directory scaffold
- [x] Add CAD/MuJoCo synchronization checks
- **Status:** complete

### Phase 5: Verification and Delivery
- [x] Run model compile and stability tests
- [x] Run ROM, scapular-rhythm, pronation, muscle-path, and flex tests
- [x] Document commands, current limitations, and next CAD import steps
- **Status:** complete

### Phase 6: Adopt Actual SolidWorks CAD Export
- [x] Locate the exported STEP/STL files
- [x] Register BodyParts3D reference bones to the exported STLs
- [x] Generate a CAD-based MuJoCo skeleton and muscle model
- [x] Revalidate flexible-body stability and muscle moment arms
- [x] Prepare the Rhino PSM generation script
- **Status:** complete

### Phase 7: Cable-driven shoulder and elbow
- [x] Add SolidWorks export macro for local STLs and named joint frames
- [x] Add metre-based mesh and joint-frame preparation
- [x] Add cable-drive manifest, fixed guides and world-fixed sliders
- [x] Add two-pass one-way tension calibration
- [x] Add compile, moment-arm, clearance and scripted-motion tests
- [ ] Run the motion acceptance test after the new five-bone STL and joint
      frames are exported
- **Status:** scaffold complete

### Phase 8: BodyParts3D arm with real muscles and printed joints
- [x] Collect the BP3D right-arm bones and derive the joint frames from geometry
- [x] Fuse 26 muscle elements into 14 channels and extract the wire centre lines
- [x] Generate the ligament strips that define each joint's range of motion
- [x] Build, check and drive `arm_bp3d.xml` and add one-click launchers
- [ ] Reach 90 degrees of shoulder flexion
- [ ] Replace the visual muscle mesh with a `flex` body and couple the tendon strain
- [ ] Concave collision decomposition so the scapulothoracic contact can be enabled
- **Status:** working model, shoulder flexion still short of target

## Key Questions
1. Exact `part.SLDASM` transforms and joint centers must be imported before final CAD validation.
2. Full-resolution flex is computationally expensive; optimization will use reduced flex grids and final candidates use higher resolution.
3. MuJoCo tendons do not collide with flex bodies; equivalent nodal-force/strain coupling is implemented.

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Right arm only, no hand | Matches the clarified user scope |
| SolidWorks is the CAD authority | Final part and assembly geometry must be manufactured |
| MuJoCo flex instead of external FEM | Keeps one Python runtime and supports rapid optimization |
| Rhino/Grasshopper for PSM generation | Matches the source paper and supports parameter gradients |
| Optuna NSGA-II for muscle search | Handles discrete muscle/subset and continuous geometry variables |
| Config-driven model generation | Prevents drift between CAD, MuJoCo, and Grasshopper |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None yet | 0 | - |
