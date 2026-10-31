# 未通过的CAD尝试

命令：`D:/Python3.10.6/python.exe mujoco_arm/anatomical_concept/export_SC_cad.py --GH --preserve-head`

实际结果：退出码1；`humerus_GH_ball`共6110个输入三角面，OpenCascade缝合/修复后仍报`RuntimeError: Invalid BRep humerus_GH_ball`。

该候选仅有网格回读及路径检查记录，不能称为已完成STEP实体交付。它还在很小的抬升角度发生实体相交，暂不采用；保留源参数及生成几何供比较。
