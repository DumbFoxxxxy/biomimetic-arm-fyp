"""Remove recorded Windows junctions non-recursively; preserve target data."""
import json
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def run():
    entries=json.loads((ROOT/'tools/workspace_links.json').read_text(encoding='utf8'));removed=[]
    for entry in entries:
        link=ROOT/entry['path'];target=(ROOT/entry['target']).resolve()
        assert link.absolute().is_relative_to(ROOT) and target.is_relative_to(ROOT)
        if not link.exists():continue
        if not link.is_junction():raise RuntimeError('Refusing ordinary directory: '+str(link))
        if link.resolve()!=target:raise RuntimeError('Unexpected junction target: '+str(link))
        os.rmdir(link)
        assert target.is_dir() and not link.exists()
        removed.append(entry)
    (ROOT/'整理记录/路径改造/removed_links.json').write_text(json.dumps(removed,ensure_ascii=False,indent=2),encoding='utf8')
    print('Removed junctions without deleting targets:',len(removed))

if __name__=='__main__':run()
