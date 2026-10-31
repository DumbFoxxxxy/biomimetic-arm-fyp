"""Copy only catalogued support references, preserving raw coordinates and sources."""
from build_reference import *
import shutil

def run():
    catalog=json.loads((OUT/'support_reference_catalog.json').read_text(encoding='utf8'));dest=ROOT/'骨骼/支撑参考';dest.mkdir(exist_ok=True);selected={}
    for category,group in catalog['structures'].items():
        for row in group['reference_models']:
            rec=selected.setdefault(row['id'],{**row,'categories':[]});rec['categories'].append(category)
    result=[]
    for fid,row in sorted(selected.items()):
        src=Path(row['source']);target=dest/src.name;sha=hashlib.sha256(src.read_bytes()).hexdigest()
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()!=sha:raise RuntimeError('Different existing user file preserved: '+str(target))
        if not target.exists():shutil.copy2(src,target)
        if hashlib.sha256(target.read_bytes()).hexdigest()!=sha:raise RuntimeError('Copy mismatch '+fid)
        result.append({**row,'workspace_source':str(target.relative_to(ROOT)),'sha256':sha,'kind':'cartilage_reference' if 'costal_cartilage' in row['categories'] else 'bone_reference','robot_part':False})
    report={'scope':'anatomical support context only; not all references must be printed','source_files_unchanged':True,'common_raw_coordinates_preserved':True,'models':result};(dest/'sources.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print('copied/verified support references',len(result),'cartilage',sum(r['kind']=='cartilage_reference' for r in result))

if __name__=='__main__':run()
