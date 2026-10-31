"""Copy verified right-side supplements into the user-designated muscle folder."""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('D:/XJTLU/isa_BP3D_4.0_obj_99')
IDS=['FJ1446','FJ1447','FJ1456','FJ1459','FJ1464','FJ1467','FJ1468','FJ1513','FJ1520','FJ1521','FJ1536','FJ1537','FJ1554']

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    dest=ROOT/'muscle';dest.mkdir(exist_ok=True);rows=[]
    for fid in IDS:
        src=SOURCE/(fid+'.obj');out=dest/src.name
        if not src.exists():raise RuntimeError('Missing source '+str(src))
        head=src.read_text(encoding='utf8',errors='replace').splitlines()[:25];name=next(s.split(' : ',1)[1] for s in head if s.startswith('# English name : '))
        if 'right' not in name.lower():raise RuntimeError('Unverified laterality '+name)
        if out.exists() and digest(out)!=digest(src):raise RuntimeError('Existing different user file preserved: '+str(out))
        copied=not out.exists()
        if copied:shutil.copy2(src,out)
        if digest(src)!=digest(out):raise RuntimeError('Copy mismatch '+fid)
        rows.append({'id':fid,'name':name,'source':str(src),'destination':str(out),'sha256':digest(out),'newly_copied':copied,'verified_identical':True})
    (dest/'supplement_sources.json').write_text(json.dumps({'source_root':str(SOURCE),'supplements':rows,'all_OBJ_count':len(list(dest.glob('*.obj'))),'source_unchanged':True},ensure_ascii=False,indent=2),encoding='utf8')
    print('supplements',len(rows),'new copies',sum(r['newly_copied'] for r in rows),'total OBJ',len(list(dest.glob('*.obj'))),'hash checks',all(r['verified_identical'] for r in rows))

if __name__=='__main__':run()
