# 50 · blache_2024_scapula_ellipsoid

- **中文标题**：椭球面是否适合建模肩胸滑动平面？
- **原题**：Is an ellipsoid surface suitable to model the scapulothoracic sliding plane?
- **作者**：Y. Blache, F. Lefebvre, I. Rogowski, B. Michaud, M. Begon
- **出处**：Journal of Biomechanics, Vol. 164, 111989（2024）
- **类型**：期刊论文　**等级**：核心　**分组**：钢丝张力分配与肩胸接触优化
- **标签**：`肩胸滑动` `单侧接触` `接触优化` `椭球接触面` `肩胛骨`
- **DOI**：10.1016/j.jbiomech.2024.111989
- **源文件**：`D:\XJTLU\PSP\literature_inbox\50_is_an_ellipsoid_surface_suitable_to_model_the_scapulothoracic_sliding_plane.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\50_is_an_ellipsoid_surface_suitable_to_model_the_scapulothoracic_sliding_plane__39e0bc.txt`
- **规模**：5 页 / 提取 28814 字符

## 摘要

论文用个体化、动作相关的椭球面拟合肩胛骨中心相对胸廓的滑动轨迹。两名受试者的平均肩胛-椭球距离分别为 0.52 mm 和 0.38 mm；只用 5 个静态姿态标定时仍为 0.57 mm 和 0.47 mm。约 5% 的姿态偏差可达 10 mm，因此作者认为硬约束并不总合适，可优先考虑允许偏离的软约束。

## 要点

- 功能标定的椭球面比单纯按胸廓外形拟合更能代表肩胸滑动平面。
- 椭球参数应针对个体和动作重新标定，不能直接套用通用球面。
- 两名受试者的平均拟合距离约 0.52 mm 和 0.38 mm。
- 仅使用 5 个静态姿态即可保持约 0.57 mm 和 0.47 mm 的中位距离。
- 局部约 5% 姿态可偏离椭球约 10 mm，硬贴面约束可能扭曲整体运动。
- 推荐使用单点或稀疏接触的软约束，而不是强制整块肩胛骨贴面。

## 与本项目的关系

当前 R400 圆柱接触带是工程简化，本文为该简化的误差评估提供了直接方法。V4 应报告接触垫自由脱离距离，并把硬穿透与允许的 1-10 mm 级偏离分开；未来可把 R400 圆柱升级为个人化椭球参数，而不必改变骨架坐标。
