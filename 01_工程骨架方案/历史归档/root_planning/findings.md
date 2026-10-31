# Findings and Decisions

## Requirements
- Right arm only.
- Model from the thorax/shoulder girdle through the distal forearm.
- No hand, finger, or wrist-DOF model.
- Use X-lattice muscle structures, not anatomical muscle meshes one-for-one.
- Use `part.SLDASM` as the eventual dimensional authority.
- SolidWorks is the main CAD environment; MuJoCo is synchronized for motion.
- Muscle count and layout are optimization variables.
- Required muscles: supraspinatus, infraspinatus, anterior/middle/posterior deltoid.
- Scapulothoracic candidates: serratus anterior, trapezius parts, rhomboids, levator scapulae, pectoralis minor.
- Elbow: two to three functional groups, including biceps/brachialis and triceps.
- Single- or dual-cable muscles and several PSM shape classes are allowed.

## Research Findings
- The Haoran Wang paper uses X-lattice TPU muscles, cable drive, 13 motor channels, a 2 mm lead screw, and fourth-order polynomial joint/muscle mappings.
- The shoulder-complex MuJoCo paper models 20 shoulder muscles and ligaments on thorax, clavicle, scapula, and humerus.
- PSM reference values are `(x,y,z,p)=(15,10,30,2) mm`; `y <= x <= z`; branch diameter strongly controls stiffness.
- The Yoshimura parametric framework uses Rhinoceros/Grasshopper and Python generation with lattice density, branch thickness, and spatial gradients.
- Wire-arrangement optimization uses a start point, relay points, an endpoint, and NSGA-II in Optuna.
- The forearm literature supports modeling pronation/supination with the radius moving around a nearly fixed ulna and with a slanted axis through the proximal radius and distal ulna.

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Build a generated MJCF rather than hand-editing one XML | Geometry, muscles, and CAD data must stay synchronized |
| Use MuJoCo `flexcomp` for PSM bodies | Keeps flexible nonlinear geometry in the main simulation |
| Use spatial tendons for cables and cable paths | Native MuJoCo tendon routing supports sites, geometry wrapping, and actuation |
| Use equivalent nodal force/strain coupling | MuJoCo tendons do not collide with flex bodies |
| Use reduced flex grids during optimization | Full lattices would make thousands of candidate evaluations impractical |
| Preserve the current pronation contact constraint | It matches the required proximal fixed / distal crossing behavior |

## Issues
| Issue | Resolution |
|-------|------------|
| Rhino/Grasshopper is not installed | Produce a documented interface and create STEP/STL placeholders |
| `part.SLDASM` is currently open | Keep current reference model and make final CAD import data-driven |
| MuJoCo environment lacked Optuna/SciPy | Installed Optuna 5.0.0 and SciPy 1.18.1 into `D:\MujocoEnv`; fallback remains available |
| Rhino path was not in `PATH` | Rhino 8 is installed at `D:\XJTLU\Rhino9.0\System`; Grasshopper is under `D:\XJTLU\Rhino9.0\Plug-ins\Grasshopper` |
| Export STLs use a transformed coordinate frame | Reference bones were registered to STLs with landmark-based Kabsch transforms |
| Hidden Rhino startup does not register RhinoCode or execute `/runscript` | Generate inside an interactive Rhino session with `_-RunPythonScript` |

## Implemented Interfaces
- `mujoco_arm/config/arm_model.json`: skeleton, ligament, gravity and CAD contract.
- `mujoco_arm/config/muscle_manifest.json`: 23-muscle candidate pool and 7 enabled required muscles.
- `mujoco_arm/build_arm_model.py`: generated MuJoCo model with flex bodies, tendons and sensors.
- `mujoco_arm/optimize_muscles.py`: Optuna NSGA-II interface with random fallback.
- `mujoco_arm/export_grasshopper_manifest.py`: JSON and CSV handoff for Grasshopper.
- `mujoco_arm/import_cad_export.py`: CAD STL registration and CAD-config generation.
- `mujoco_arm/grasshopper_interface/generate_psm_rhino.py`: interactive Rhino PSM generation and export.

## CAD Registration Results

The following files were generated from the SolidWorks STL exports:

- `mujoco_arm/arm_cad_reference.xml`
- `mujoco_arm/arm_psm_generated_cad.xml`
- `mujoco_arm/config/arm_model_cad.json`
- `mujoco_arm/config/muscle_manifest_cad.json`
- `mujoco_arm/config/cad_registration.json`

Reference-to-export sample RMSE:

| Bone | RMSE |
|------|------|
| clavicle | 5.682 mm |
| scapula | 7.943 mm |
| humerus | 12.511 mm |
| radius | 5.919 mm |
| ulna | 0.577 mm |

## Resources
- `D:\XJTLU\PSP\2141387_Haoran Wang_2025.pdf`
- `D:\XJTLU\PSP\Patterned Structure Muscle  Arbitrary Shaped Wire-driven Artificial Muscle Utilizing Anisotropic Flexible Structure for Musculoskeletal Robots.pdf`
- `D:\XJTLU\PSP\Advanced Intelligent Systems - 2026 - Yoshimura - Accelerating Musculoskeletal Robotics Through Parametric Design and.pdf`
- `D:\XJTLU\PSP\Construction_of_Musculoskeletal_Simulation_for_Shoulder_Complex_with_Ligaments_and_Its_Validation_via_Model_Predictive_Control.pdf`
- `D:\XJTLU\PSP\Design_Optimization_of_Wire_Arrangement_With_Variable_Relay_Points_in_Numerical_Simulation_for_Tendon-Driven_Robots.pdf`
- `D:\XJTLU\PSP\mujoco_arm`

## Visual Findings
- PSM consists of an X lattice, rigid wire-channel blocks, an outer shell, and dense end connectors.
- The shoulder includes thorax, clavicle, scapula, and humerus; scapulothoracic motion requires thorax/rib attachment regions.
- The elbow region needs humeroulnar, humeroradial, proximal radioulnar, and distal radioulnar behavior.
