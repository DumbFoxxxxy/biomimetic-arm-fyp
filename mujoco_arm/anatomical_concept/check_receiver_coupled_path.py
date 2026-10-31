"""Cover actual combined yaw/pitch poses, beyond unit radial samples."""
from build_forearm_receivers import *
MOUNT=DEST/'rephased_anchors'

def run():
    source=MOUNT/'mounted_direction_screen.json';screen=json.loads(source.read_text());reports=[]
    for fid in ['FJ1486','FJ1477']:
        d=json.loads((DEST/fid/'design.json').read_text());parts={k:load(DEST/fid/(k+'.obj'))[0] for k in d['parts']};first=None
        for row in screen['poses']:
            angles=row['modules'][fid];posed=articulate(parts,angles['yaw_relative_lever_deg'],angles['pitch_deg'],d['standoff_mm']);h=[]
            for (a,sa),(b,sb) in itertools.combinations(posed.items(),2):
                v=intersection(sa,sb)
                if v>.001:h.append({'pair':[a,b],'volume_mm3':v})
            if h and first is None:first={'elbow_angle_deg':row['elbow_angle_deg'],'hits':h}
        reports.append({'muscle':fid,'coupled_samples':len(screen['poses']),'first_internal_limit':first,'source_sha256':{k:hashlib.sha256((DEST/fid/(k+'.obj')).read_bytes()).hexdigest() for k in parts}})
    result={'scope':'actual combined local yaw/pitch configurations from conditional mounted direction study, all13-part pairs; not a valid complete tendon route','modules':reports,'all_clear':not any(r['first_internal_limit'] for r in reports),'screen_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'frozen':False};(MOUNT/'coupled_path_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M8 coupled path clear',result['all_clear'],flush=True);assert result['all_clear']

if __name__=='__main__':run()
