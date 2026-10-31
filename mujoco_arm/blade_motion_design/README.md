# Blade paper shoulder motion design

> Scope correction, 2026-10-29: the user requires wire + TPU patterned artificial muscles. The current models contain ideal tendons and rigid geometry only. Their passing motion checks do not validate TPU muscle actuation and do not complete the FYP drive design. See `../../design/TPU_MUSCLE_SCOPE_CORRECTION.md` before reusing anchors, cable counts or the concentric GH candidate.

This directory implements the active shoulder-design objective, not the earlier proposed five-active-DOF arm simplification. It retains SC2 + AC3 and tests shoulder-blade motion with eight pull-only tendons. Payload and structural strength are not acceptance criteria. No native CAD file is overwritten.

The main paper is `A_shoulder_structure_of_muscle-driven_humanoid_with_shoulder_blades.pdf`. Four crossed upper/lower routes follow its page-5 topology; four additional routes stabilize this project's five-DOF shoulder. The final coordinates are derived for this project, not copied from the paper.

Current design uses a 12 mm outward scapula-body stand-off and an 8 mm anterior clavicle-body offset, with webs reconnecting the original joint centers. A genuine boolean-subtracted S07 corridor permits compact local routing. Other cables use exterior eyes. The final 3 mm at each cable termination is a designated clamp/eye contact region; it is not a freely sliding channel.

## Current scope

`../exports/blade_motion_design/shoulder_checkpoint.json` and `shoulder_anchor_checkpoint.csv` are a shoulder-subsystem checkpoint. They are not frozen full-arm manufacturing hardpoints. Full GH routing, full-arm interactions and final ROM remain in progress.

## Run the shoulder demo

From the project root:

```powershell
& 'D:\MujocoEnv\Scripts\python.exe' mujoco_arm/blade_motion_design/run_shoulder_demo.py
```

Keys 1-6 queue up/down rotation, protraction/retraction and elevation/depression for the next cycle. No state reset occurs after initialization. Motion is driven by eight tendon actuators, not prescribed joint angles. The current diagnostic model uses no gravity and no payload.

For a short non-GUI run:

```powershell
& 'D:\MujocoEnv\Scripts\python.exe' mujoco_arm/blade_motion_design/run_shoulder_demo.py --headless --seconds 1
```

## Verification

- `audit_shoulder_drive.py`: actual tendon-driven tracking, including return motion.
- `audit_geometry.py`: dense body/chest geometry and wire centerline/mesh checks.
- `audit_channel_clearance.py`: finite wire radius in the S07 free corridor and pairwise wire separation, using independent segment/triangle distances.
- `audit_arm_rom.py`: integrated GH/elbow geometric exploration only; it does not prove GH tendon control or all self-collisions.
- `export_shoulder_checkpoint.py`: verifies input fingerprints match all shoulder reports before exporting coordinates.

The all-compact eight-channel option is retained as a rejected experiment: large swept openings split the plate into multiple substantial pieces. It is not substituted into the accepted local-channel model.

Boolean construction uses manifold3d installed into the project-local `tmp/shoulder_motion_design_20261029/vendor`. Existing V8 code supplies coordinate/mesh helpers. Exported model and reports include absolute local paths or input hashes for reproducibility.
