"""Anatomical structure constraints, independent from old cable/ROM designs.
Atlas identifies structures; geometric patches are separate review candidates.
"""
UW='https://www.rad.washington.edu/muscle-atlas/'
BONES={'clavicle':'FJ3362','scapula':'FJ3384','humerus':'FJ3368','ulna':'FJ3391','radius':'FJ3349'}

def atlas():
    rows={}
    def add(ids,origin,insertion,origin_region,insertion_region,source,broad=False):
        for fid in ids.split():rows[fid]={'origin_structures':origin.split(),'insertion_structures':insertion.split(),'origin_named_region':origin_region,'insertion_named_region':insertion_region,'source_url':source if source.startswith('http') else UW+source,'source_review_date':'2026-10-30','broad_attachment_review':broad,'anatomical_landmark_coordinates_verified':False,'drill_or_wire_anchor_frozen':False}
    add('FJ1460','first_rib costal_cartilage','clavicle','第一肋与肋软骨交界','锁骨中段下方沟','https://www.meddean.luc.edu/lumen/meded/grossanatomy/dissector/mml/subc.htm')
    add('FJ1477','humerus','ulna','肱骨后侧、桡神经沟上方','尺骨鹰嘴近端','triceps-brachii',True)
    add('FJ1479','scapula','ulna','肩胛盂下结节','尺骨鹰嘴近端','triceps-brachii')
    add('FJ1480','humerus','ulna','肱骨后侧、桡神经沟下方','尺骨鹰嘴近端','triceps-brachii',True)
    add('FJ1478','scapula','radius','肩胛盂上结节','桡骨粗隆；另有腱膜连接不在本骨模型中','biceps-brachii')
    add('FJ1512','scapula','radius','肩胛喙突尖部','桡骨粗隆；另有腱膜连接不在本骨模型中','biceps-brachii')
    add('FJ1486','humerus','ulna','肱骨远端前侧','尺骨冠突与粗隆','brachialis',True)
    add('FJ1488','scapula','humerus','肩胛喙突尖部','肱骨内侧中段','coracobrachialis')
    add('FJ1500','scapula','humerus','肩胛冈下窝','肱骨大结节中面','infraspinatus',True)
    add('FJ1504','scapula','humerus','肩胛下窝','肱骨小结节','subscapularis',True)
    add('FJ1506','scapula','humerus','肩胛冈上窝','肱骨大结节上面','supraspinatus',True)
    add('FJ1507','scapula','humerus','肩胛下角后侧','肱骨结节间沟内侧唇','teres-major')
    add('FJ1508','scapula','humerus','肩胛外侧缘上部','肱骨大结节下面','teres-minor')
    add('FJ1532','cervical_spine','scapula','C1–C4横突区域','肩胛内侧缘上段','levator-scapulae')
    add('FJ1468','clavicle','humerus','锁骨外侧段','肱骨三角肌粗隆','deltoid',True)
    add('FJ1467','scapula','humerus','肩峰','肱骨三角肌粗隆','deltoid',True)
    add('FJ1513','scapula','humerus','肩胛冈','肱骨三角肌粗隆','deltoid',True)
    add('FJ1447','clavicle','humerus','锁骨前侧内段','肱骨结节间沟外侧唇','pectoralis-major',True)
    add('FJ1464','sternum costal_cartilage ribs','humerus','胸骨、上部肋软骨相关区域','肱骨结节间沟外侧唇','pectoralis-major',True)
    add('FJ1446','abdominal_aponeurosis','humerus','胸大肌腹部的腱膜起始区域','肱骨结节间沟外侧唇','pectoralis-major',True)
    add('FJ1456','ribs costal_cartilage','scapula','第三至第五肋附近','喙突上内侧','pectoralis-minor',True)
    add('FJ1459','ribs','scapula','胸侧上八至九肋相关区域','肩胛内侧缘','serratus-anterior',True)
    add('FJ1521 FJ1554 FJ1520','occiput nuchal_ligament cervical_spine thoracic_spine','clavicle scapula','斜方肌整体起始区；各分部区域须按模型另行标定','斜方肌整体止区含锁骨外段、肩峰和肩胛冈；不要求每分部同时连接两骨','trapezius',True)
    add('FJ1536','thoracic_spine','scapula','T2–T5棘突区域','肩胛内侧缘相关区域','rhomboid-major-and-minor',True)
    add('FJ1537','nuchal_ligament cervical_spine thoracic_spine','scapula','项韧带及C7、T1区域','肩胛内侧缘相关区域','rhomboid-major-and-minor')
    return rows
