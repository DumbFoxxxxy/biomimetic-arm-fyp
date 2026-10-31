"""Limited straight pullout probes for the AC candidate, not a retention proof."""
from build_AC_connection import *

def run():
    data=json.loads((DEST/'assembly.json').read_text(encoding='utf8'))
    centre=np.array(data['centre_source_mm'])-SHIFT
    z=np.array(data['socket_axis_toward_neck']);x=unit(np.cross(z,[0,0,1]));y=np.cross(z,x)
    fixed={}
    for name in ['scapula_AC_socket','AC_cap_left','AC_cap_right']:
        v,f=read_obj(str(DEST/(name+'.obj')))
        fixed[name]=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
    ball=mf.Manifold.sphere(data['ball_radius_mm'],32).translate(centre);rows=[]
    for name,axis in [('opening',z),('rear',-z),('side_x',x),('side_minus_x',-x),('side_y',y),('side_minus_y',-y)]:
        first=None
        for distance in np.arange(.1,2.001,.1):
            moved=ball.translate(axis*distance)
            hits=[{'part':n,'overlap_mm3':float((moved^s).volume())} for n,s in fixed.items()]
            hits=[h for h in hits if h['overlap_mm3']>.001]
            if hits:first={'translation_mm':float(distance),'contacts':hits};break
        rows.append({'direction':name,'first_sampled_interference':first})
    report={'scope':'ball-only straight translations; socket/caps fixed. Not arbitrary-path escape proof or retention-load test.','step_mm':.1,'maximum_probe_mm':2,'directions':rows,'all_six_straight_directions_blocked_within_2mm':all(r['first_sampled_interference'] is not None for r in rows),'requires':'fasteners, printed clearance and strength validation'}
    (DEST/'capture_probes.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('AC straight pullout probes',[(r['direction'],r['first_sampled_interference']['translation_mm'] if r['first_sampled_interference'] else None) for r in rows])

if __name__=='__main__':run()
