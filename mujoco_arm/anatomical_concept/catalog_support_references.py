"""Locate anatomical support references without treating them as robot fixtures."""
from build_reference import *
import csv

def run():
    base=Path('D:/XJTLU/isa_BP3D_4.0_obj_99');index=base/'isa_element_parts.txt';rows=list(csv.DictReader(index.open(encoding='utf8'),delimiter='\t'));candidates={}
    ordinals='first|second|third|fourth|fifth|sixth|seventh|eighth|ninth'
    patterns={'first_rib':r'^right first rib$','ribs':rf'^right ({ordinals}) rib$','costal_cartilage':rf'^right ({ordinals}) costal cartilage$','sternum':r'^(manubrium|body of sternum|xiphoid process)$','cervical_spine':r'^(atlas|axis|(first|second|third|fourth|fifth|sixth|seventh) cervical vertebra)$','thoracic_spine':r'^(first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|eleventh|twelfth) thoracic vertebra$','occiput':r'^occipital bone$','nuchal_ligament':r'^nuchal ligament$','abdominal_aponeurosis':r'^(aponeurosis of external abdominal oblique|external oblique aponeurosis)$'}
    for category,pattern in patterns.items():
        hits=[];seen=set()
        for r in rows:
            fid=r['element file id'];p=base/(fid+'.obj')
            if re.search(pattern,r['name'],re.I) and fid not in seen and not fid.endswith('M') and p.exists():hits.append({'id':fid,'label':r['name'],'concept':r['concept id'],'source':str(p)});seen.add(fid)
        candidates[category]={'reference_models':hits,'status':'references available for later mechanical registration' if hits else 'no matching model label identified; soft-tissue or frame reference needs explicit definition','robot_fixture_designed':False}
    report={'scope':'available source reference geometry, not installed support frames','source_index':str(index),'index_sha256':hashlib.sha256(index.read_bytes()).hexdigest(),'structures':candidates};(OUT/'support_reference_catalog.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print({k:len(v['reference_models']) for k,v in candidates.items()})

if __name__=='__main__':run()
