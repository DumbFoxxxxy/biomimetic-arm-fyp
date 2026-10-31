# TPU / steel-wire muscle design — active work

> Latest user direction: skeleton validation from V2 comes first. Coupon work is retained as a later interface-design reference; do not continue coupon optimization ahead of the skeletal baseline. See `../../design/V2_SKELETON_FIRST_SCOPE.md`.

This replaces bare-tendon actuation as the final FYP architecture. Original CAD and earlier ideal-cable studies are retained as references only.

## Current artifact

`../exports/tpu_muscle_design/psm_reference_coupon_mm.stl` is an **undeformed internal-structure coupon**, in millimetres. It uses a Figure-5-inspired centre-to-eight-corners abstraction, not the author's original CAD. Five 15 × 10 × 30 mm cells have 2 mm branches. Discrete wire-guide rings and 10 mm dense end blocks are connected. The overall envelope is 95 × 12 × 32 mm; the 75 mm cell span excludes the two end blocks. The central guide bore is 1.8 mm. End fixtures have two 3.2 mm clearance holes each.

The manifold build verifies one material component and a clear axial wire passage. The preview was inspected. These checks do **not** establish printable overhang quality, contraction, strain safety, cyclic life, force output, or compatibility with the arm. The outer patterned surface is not implemented yet.

Run with the project's existing dependencies:

```powershell
& 'D:\MujocoEnv\Scripts\python.exe' mujoco_arm/tpu_muscle_design/build_reference_coupon.py
& 'C:\Python314\python.exe' mujoco_arm/tpu_muscle_design/render_coupon.py
```

## Stroke evidence and assumptions

`reference_stroke_budget.csv` screens the earlier cable candidate; its anchor coordinates are not adopted for TPU. The columns Δwire / assumed contraction fraction apply only to a hypothetical 1:1 shortening transmission, exclude fixed ends, and do not certify material behavior. They intentionally show scenarios at 18%, 20%, 35%, and 50%, rather than inventing a validated strain limit.

For straight shoulder routes there is also an illustrative rest-length interval based on the PSM paper's **whole muscle** length ratios: nominal 250 mm, contracted 115 mm, stretched 380 mm (Table I). This accounts for extension as well as shortening: Lrest ≥ Lmax / 1.52 and Lrest ≤ Lmin / 0.46. These are literature-scaled screening bounds, not limits measured for this coupon. They cannot be applied to the GH/elbow cable lengths, which include wrap arcs. A fixed-geometry wire stroke is not automatically the same as TPU shortening.

The AIS2026 paragraph around Figure 4 calls ~20% a biological-muscle macroscopic reference. It does not impose a universal 20% PSM ceiling. The old design brief's attribution is corrected accordingly.

## Required next evidence

### Compression-assay progress

`solve_coupon_rods.py` now minimizes stretching, discrete bending and junction-angle energies with end displacement imposed only at the fixtures. Internal nodes are solved, not affinely scaled. The eight-corner branch topology has 149 nodes and 160 rod segments; straight wire guides constrain central-node lateral motion. Young's modulus is factored out consistently, so the model predicts a normalized equilibrium geometry, not forces in newtons. An independent directional gradient check has absolute error below 6e-11. All 21 sampled states satisfy the declared free-gradient residual threshold 1e-6; optimizer status alone is not used as proof of equilibrium.

At 25% cell-span shortening, the span changes 75 -> 56.25 mm. With the two rigid 10 mm end blocks, this is 95 -> 76.25 mm overall (19.74% overall shortening). The rod-centreline transverse envelope grows from 10 × 30 to approximately 11.71 × 31.31 mm, before branch radii and any outer layer are added. This change in occupied volume must be included in shoulder installation checks.

The first sampled new nonjunction branch intersection occurs at 30% cell-span shortening. The contact screen trims 4 mm around original fused junctions and therefore is not a full self-contact certificate. The solve has no contact forces: states after detected contact are invalid as predictions of the real material. The 30% observation is neither an experimental contraction limit nor a validated operating rating. Outer-pattern constraints, wire-guide-ring stiffness, junction torsion, printing anisotropy, large-strain material calibration and equilibrium stability remain unresolved. The current result is a compression assay, not evidence that the muscle drives a skeletal joint.

Outputs: `coupon_rod_compression_report.json`, `coupon_rod_equilibria.npz`, and inspected `coupon_rod_assay.png`.

```powershell
& 'D:\MujocoEnv\Scripts\python.exe' mujoco_arm/tpu_muscle_design/solve_coupon_rods.py
& 'C:\Python314\python.exe' mujoco_arm/tpu_muscle_design/render_rod_assay.py
```

1. Complete the guide and outer-pattern geometry and evaluate compression, transverse deformation, branch self-contact and wire sliding. An affine mesh scaling or prescribed animation is insufficient proof of muscle-driven motion.
2. Establish the actual shortening / extension envelope using a mechanically meaningful model and/or a physical coupon; no measured filament constitutive data currently exists in this workspace.
3. Couple an antagonistic muscle pair to the elbow, with both muscle endpoints attached to the skeletal links and the internal wire included.
4. Fit muscle shape and installation to the shoulder, preserving the user's full motion targets. Only then derive final anchors and ROM. The direct-cable layout, 80 mm GH relocation, and 17-line count are not frozen constraints.

No payload/strength acceptance is required. Contraction kinematics and occupied TPU volume remain required.
