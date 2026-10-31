# 分类依据

方案一：原工程骨架、V2–V8、V2-S1、旧钢丝/PSM路线；final_model、model、既有archive属于这一历史链。
方案二：真实人体骨骼与异形TPU路线；anatomical_arm为早期探索，anatomical_concept为当前肘部候选。
骨骼、muscle是方案二输入。两条方案不共享坐标、锚点或通过结论。
knowledge_base、literature_inbox、tmp含共用检索与依赖；保留现有运行位置，建立共用入口。
design/DESIGN_BRIEF.md混合历史与当前记录，保留为历史总账，单独建立两条路线的文档入口。
脚本大量使用Path(__file__).resolve().parents[2]，移动源代码会改变根目录推导，因此只迁移数据，代码用目录联接按路线归类。
根目录旧PROJECT_STRUCTURE和HANDOFF偏向旧路线，需保留旧稿并更新导航，避免误读。
