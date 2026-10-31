"""Finite centre-translation probes of the spherical GH head against real cage."""
from build_GH_connection import *

def run():
    data=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));G=np.array(data['centre_source_mm'])-SHIFT
    fixed={}
    for name in ['scapula_AC_GH_sockets',*[n for n in data['components'] if n.startswith('GH_cap_')]]:
        fixed[name]=load(DEST/(name+'.obj'))[0]
    ball=mf.Manifold.sphere(data['ball_radius_mm'],40).translate(G);rows=[]
    for rawaxis in itertools.product([-1,0,1],repeat=3):
        if not any(rawaxis):continue
        axis=unit(rawaxis);first=None
        for distance in np.arange(.1,3.001,.1):
            moved=ball.translate(axis*distance)
            hits=[{'part':name,'volume_mm3':float((moved^s).volume())} for name,s in fixed.items()]
            hits=[r for r in hits if r['volume_mm3']>.001]
            if hits:first={'translation_mm':float(distance),'hits':hits};break
        rows.append({'unit_direction':axis.tolist(),'first_interference':first})
    report={'scope':'26 straight head-centre translations with cage/caps held fixed; not arbitrary escape-path proof, neck assembly proof or retention-load test','head_radius_mm':data['ball_radius_mm'],'step_mm':.1,'maximum_probe_mm':3,'directions':rows,'all_26_directions_blocked_within_3mm':all(r['first_interference'] is not None for r in rows),'frozen':False}
    (DEST/'capture_probes.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('GH straight pullouts blocked',sum(r['first_interference'] is not None for r in rows),'of',len(rows),'max sampled movement',max((r['first_interference']['translation_mm'] for r in rows if r['first_interference']),default=None))

if __name__=='__main__':run()
