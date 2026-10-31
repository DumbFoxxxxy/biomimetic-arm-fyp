# 肌肉附着关系与几何候选区域

以下结构对应关系按大学解剖图谱整理。几何区域只是在这些指定骨骼上做表面距离筛查，
不等同于已确认的解剖附着区、钢丝锚点或可钻孔坐标。PCA方向没有用于推断肌纤维。

|ID|肌肉|起始结构|止端结构|参考来源|
|---|---|---|---|---|
|FJ1446|Abdominal part of right pectoralis major|胸大肌腹部的腱膜起始区域|肱骨结节间沟外侧唇|[atlas](https://www.rad.washington.edu/muscle-atlas/pectoralis-major)|
|FJ1447|Clavicular part of right pectoralis major|锁骨前侧内段|肱骨结节间沟外侧唇|[atlas](https://www.rad.washington.edu/muscle-atlas/pectoralis-major)|
|FJ1456|Right pectoralis minor|第三至第五肋附近|喙突上内侧|[atlas](https://www.rad.washington.edu/muscle-atlas/pectoralis-minor)|
|FJ1459|Right serratus anterior|胸侧上八至九肋相关区域|肩胛内侧缘|[atlas](https://www.rad.washington.edu/muscle-atlas/serratus-anterior)|
|FJ1460|Right subclavius|第一肋与肋软骨交界|锁骨中段下方沟|[atlas](https://www.meddean.luc.edu/lumen/meded/grossanatomy/dissector/mml/subc.htm)|
|FJ1464|Sternocostal part of right pectoralis major|胸骨、上部肋软骨相关区域|肱骨结节间沟外侧唇|[atlas](https://www.rad.washington.edu/muscle-atlas/pectoralis-major)|
|FJ1467|Acromial part of right deltoid|肩峰|肱骨三角肌粗隆|[atlas](https://www.rad.washington.edu/muscle-atlas/deltoid)|
|FJ1468|Clavicular part of right deltoid|锁骨外侧段|肱骨三角肌粗隆|[atlas](https://www.rad.washington.edu/muscle-atlas/deltoid)|
|FJ1477|Lateral head of right triceps brachii|肱骨后侧、桡神经沟上方|尺骨鹰嘴近端|[atlas](https://www.rad.washington.edu/muscle-atlas/triceps-brachii)|
|FJ1478|Long head of right biceps brachii|肩胛盂上结节|桡骨粗隆；另有腱膜连接不在本骨模型中|[atlas](https://www.rad.washington.edu/muscle-atlas/biceps-brachii)|
|FJ1479|Long head of right triceps brachii|肩胛盂下结节|尺骨鹰嘴近端|[atlas](https://www.rad.washington.edu/muscle-atlas/triceps-brachii)|
|FJ1480|Medial head of right triceps brachii|肱骨后侧、桡神经沟下方|尺骨鹰嘴近端|[atlas](https://www.rad.washington.edu/muscle-atlas/triceps-brachii)|
|FJ1486|Right brachialis|肱骨远端前侧|尺骨冠突与粗隆|[atlas](https://www.rad.washington.edu/muscle-atlas/brachialis)|
|FJ1488|Right coracobrachialis|肩胛喙突尖部|肱骨内侧中段|[atlas](https://www.rad.washington.edu/muscle-atlas/coracobrachialis)|
|FJ1500|Right infraspinatus muscle|肩胛冈下窝|肱骨大结节中面|[atlas](https://www.rad.washington.edu/muscle-atlas/infraspinatus)|
|FJ1504|Right subscapularis|肩胛下窝|肱骨小结节|[atlas](https://www.rad.washington.edu/muscle-atlas/subscapularis)|
|FJ1506|Right supraspinatus|肩胛冈上窝|肱骨大结节上面|[atlas](https://www.rad.washington.edu/muscle-atlas/supraspinatus)|
|FJ1507|Right teres major|肩胛下角后侧|肱骨结节间沟内侧唇|[atlas](https://www.rad.washington.edu/muscle-atlas/teres-major)|
|FJ1508|Right teres minor|肩胛外侧缘上部|肱骨大结节下面|[atlas](https://www.rad.washington.edu/muscle-atlas/teres-minor)|
|FJ1512|Short head of right biceps brachii|肩胛喙突尖部|桡骨粗隆；另有腱膜连接不在本骨模型中|[atlas](https://www.rad.washington.edu/muscle-atlas/biceps-brachii)|
|FJ1513|Spinal part of right deltoid|肩胛冈|肱骨三角肌粗隆|[atlas](https://www.rad.washington.edu/muscle-atlas/deltoid)|
|FJ1520|Ascending part of right trapezius|斜方肌整体起始区；各分部区域须按模型另行标定|斜方肌整体止区含锁骨外段、肩峰和肩胛冈；不要求每分部同时连接两骨|[atlas](https://www.rad.washington.edu/muscle-atlas/trapezius)|
|FJ1521|Descending part of right trapezius|斜方肌整体起始区；各分部区域须按模型另行标定|斜方肌整体止区含锁骨外段、肩峰和肩胛冈；不要求每分部同时连接两骨|[atlas](https://www.rad.washington.edu/muscle-atlas/trapezius)|
|FJ1532|Right levator scapulae|C1–C4横突区域|肩胛内侧缘上段|[atlas](https://www.rad.washington.edu/muscle-atlas/levator-scapulae)|
|FJ1536|Right rhomboid major|T2–T5棘突区域|肩胛内侧缘相关区域|[atlas](https://www.rad.washington.edu/muscle-atlas/rhomboid-major-and-minor)|
|FJ1537|Right rhomboid minor|项韧带及C7、T1区域|肩胛内侧缘相关区域|[atlas](https://www.rad.washington.edu/muscle-atlas/rhomboid-major-and-minor)|
|FJ1554|Transverse part of right trapezius|斜方肌整体起始区；各分部区域须按模型另行标定|斜方肌整体止区含锁骨外段、肩峰和肩胛冈；不要求每分部同时连接两骨|[atlas](https://www.rad.washington.edu/muscle-atlas/trapezius)|

## 尚缺的支撑定义

当前仅有五块上肢骨骼。以下结构并没有被替换成虚构固定点：
abdominal_aponeurosis, nuchal_ligament

后续应为这些区域建立真实参考或明确的工程机架连接。
斜方肌分部的整体附着范围仅用于候选筛查，不要求每个分部同时连接锁骨与肩胛骨。
原点和止端的骨骼对应关系确定后，还需在模型上定位具体骨性标志，并设计分布式锚定。

NPZ保留所有被筛选的肌肉面、骨面投影点和距离。标记为绿色起始区域、橙色止端区域。
显示20个抽样点不意味着只有20个锚点，也不意味着所有候选点均可制造。
