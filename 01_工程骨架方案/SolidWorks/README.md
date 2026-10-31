# Final CAD files

2026-10-29：新选定的 **V2-S1骨架几何基线**位于 `../mujoco_arm/exports/skeleton_from_v2/V2-S1/`，详见 `../design/V2_SKELETON_FREEZE.md`。新交付是经回读复核的可编辑分面BRep STEP；本目录原生SolidWorks零件/装配未回写，不能视为已与V2-S1同步。

当前保留：

```text
Thorax.SLDPRT
clavicle.SLDPRT
scapula.SLDPRT
humerus.SLDPRT
ulna.SLDPRT
universal_joint_pin.SLDPRT
assembly_zero.SLDASM
assembly_motion_gh.SLDASM
```

- `assembly_zero.SLDASM`：当前干净零位装配。
- `assembly_motion_gh.SLDASM`：最新 GH Motion 实验版本。

历史装配、旧的 14 路径模型和一次性宏在 `archive/`。
