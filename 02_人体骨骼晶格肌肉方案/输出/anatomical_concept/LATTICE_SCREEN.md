# 解剖体积晶格筛查

两组结构使用同一人体肌肉体积、同一10×10×15mm胞元和2mm名义杆径。
PCA仅用于铺设网格坐标，不被解释为肌纤维方向。第二组增加由源表面简化得到的
三角形外网格，所有材料均裁切在原体积内。这是几何筛查，不是论文曲面参数化算法
的完整复现，也没有钢丝通道或安装端。

|ID|肌肉|源体积体块|仅BCC体块|加外网格体块|材料体积分数|
|---|---|---:|---:|---:|
|FJ1460|Right subclavius|2|14|2|56.2%|
|FJ1477|Lateral head of right triceps brachii|1|11|1|23.2%|
|FJ1478|Long head of right biceps brachii|1|25|1|23.7%|
|FJ1479|Long head of right triceps brachii|1|10|1|22.4%|
|FJ1480|Medial head of right triceps brachii|1|12|1|25.0%|
|FJ1486|Right brachialis|1|21|1|27.6%|
|FJ1488|Right coracobrachialis|1|4|1|29.6%|
|FJ1500|Right infraspinatus muscle|1|10|1|20.1%|
|FJ1504|Right subscapularis|1|25|1|29.4%|
|FJ1506|Right supraspinatus|2|14|2|25.8%|
|FJ1507|Right teres major|2|8|2|26.1%|
|FJ1508|Right teres minor|1|8|1|24.7%|
|FJ1512|Short head of right biceps brachii|1|17|1|26.2%|
|FJ1532|Right levator scapulae|1|18|1|54.7%|

3个源肌肉本身带有独立的微小闭合体：它们不是晶格新产生的断裂。
这些体块在本轮保留并报告；不能为了单一体块通过而悄悄删除。
只考察体块数量仍不足以证明打印可行，须结合最细杆宽、纤维走向、
钢丝孔道、端部过渡和实际变形。

OBJ以双精度坐标保存，STL采用浮点格式，分别回读检查；状态见JSON。
STL不通过的模型仅作为查看预览，不能作为制造通过的证据。
材料体积分数较高的薄片肌肉可能过硬，不能把连通性改善当成柔性性能已验证。
