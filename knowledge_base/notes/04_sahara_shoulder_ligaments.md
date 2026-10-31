# 04 · sahara_shoulder_ligaments

- **中文标题**：含韧带的肩复合体肌肉骨骼仿真构建及其基于模型预测控制的验证
- **原题**：Construction of Musculoskeletal Simulation for Shoulder Complex with Ligaments and Its Validation via Model Predictive Control
- **作者**：Yuta Sahara, Akihiro Miki, Yoshimoto Ribayashi, Shunnosuke Yoshimura, Kento Kawaharazuka, Kei Okada, Masayuki Inaba
- **出处**：IEEE 会议论文（PDF 未标注具体会议）（2021）
- **类型**：会议论文　**等级**：核心　**分组**：项目直接依据
- **标签**：`肩复合体` `韧带` `锁骨` `肩胛骨` `MPC` `最大肌力` `仿真`
- **源文件**：`D:\XJTLU\PSP\Construction_of_Musculoskeletal_Simulation_for_Shoulder_Complex_with_Ligaments_and_Its_Validation_via_Model_Predictive_Control.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Construction_of_Musculoskeletal_Simulation_for_Shoulder_Complex_with_Ligaments_and_Its_Val__53ddd0.txt`
- **规模**：7 页 / 提取 32223 字符

## 摘要

指出现有肌肉骨骼仿真器只建肌肉不建韧带，也没有专门针对包含锁骨与肩胛骨的肩复合体。本文构建了含韧带的肩复合体详细仿真模型，同时复现骨骼结构、肌肉排列、韧带排列和最大肌力，再用模型预测控制验证。结论是韧带有助于第一个动作中关节的稳定，而最大肌力的合理分配有助于均衡各肌肉负载。

## 要点

- 明确指出「只建肌肉不建韧带」是现有仿真器的缺陷，这是加入 81 条韧带条带的直接理由。
- 肩复合体必须包含锁骨和肩胛骨，不能只做盂肱关节。
- 每条肌肉的最大力是需要标定的设计量，分配是否合理直接决定负载均衡。
- 用 MPC 做验证：给定动作反推需要多少肌肉张力。

## 与本项目的关系

你在早期方案里把它列为肩部肌肉与韧带分组的起点。韧带 manifest 与关节刚度反推思路与之一致。
