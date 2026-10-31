# 01 · psm_yoshimura_iros2020

- **中文标题**：图案化结构肌肉（PSM）：利用各向异性柔性结构的任意形状钢丝驱动机器人肌肉
- **原题**：Patterned Structure Muscle: Arbitrary Shaped Wire-driven Artificial Muscle Utilizing Anisotropic Flexible Structure for Musculoskeletal Robots
- **作者**：Shunnosuke Yoshimura, Akihiro Miki, Kazuhiro Miyama, Yuta Sahara, Kento Kawaharazuka, Kei Okada, Masayuki Inaba
- **出处**：IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)（2020）
- **类型**：会议论文　**等级**：核心　**分组**：项目直接依据
- **标签**：`X晶格` `PSM` `TPU` `FDM打印` `钢丝驱动` `各向异性`
- **源文件**：`D:\XJTLU\PSP\Patterned Structure Muscle  Arbitrary Shaped Wire-driven Artificial Muscle Utilizing Anisotropic Flexible Structure for Musculoskeletal Robots.pdf`
- **提取全文**：`D:\XJTLU\PSP\knowledge_base\fulltext\Patterned_Structure_Muscle_Arbitrary_Shaped_Wire_driven_Artificial_Muscle_Utilizing_Anisot__859f01.txt`
- **规模**：8 页 / 提取 39938 字符

## 摘要

PSM 就是你所说的 X 晶格仿生肌肉的原始出处。做法是用 FDM 打印的 TPU 图案化结构产生各向异性：沿长度方向抗压、沿横向可弯，内部穿一根钢丝；钢丝收紧时结构被压缩而变短，模拟肌肉收缩。论文给出四类形态：单自由度肌肉、多自由度宽幅肌肉、包覆关节的肌肉、Y 形分支肌肉，并用一条上臂验证了大幅运动、举重与环境接触中的运动。

## 要点

- 肌肉本体是 TPU 打印的各向异性图案结构，钢丝只负责传力，收缩变形由结构本身决定。
- 形态可自由定制，这正是你需要确定肌肉数量、位置、形状和晶格密度的原因。
- 含上臂结构实验：证明该方案能完成大范围运动、重物举升与环境接触。
- X 晶格的密度与尺寸直接决定刚度与行程，是后续要优化的设计变量。

## 与本项目的关系

这是你方案的第一原理来源。当前 arm_bp3d 管线里钢丝走肌肉体积中心线，正是为了实现这篇的实物布置：钢丝在 X 晶格内部。
