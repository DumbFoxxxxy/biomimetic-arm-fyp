"""Six straight pullout probes with the socket/caps rigidly held in place."""
from build_SC_connection import *

def run():
    definition=json.loads((BASE/'definition.json').read_text(encoding='utf8'));S=np.array(definition['centres_source_mm']['SC'])-SHIFT;A=np.array(definition['centres_source_mm']['AC'])-SHIFT;z=unit(A-S);x=unit(np.cross(z,[0,0,1]));y=np.cross(z,x);fixed={}
    for name in ['manubrium_socket_base','SC_cap_left','SC_cap_right']:
        v,f=read_obj(str(DEST/(name+'.obj')));fixed[name]=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
    ball=mf.Manifold.sphere(10,32).translate(S);rows=[]
    for name,axis in [('opening',z),('rear',-z),('side_x',x),('side_minus_x',-x),('side_y',y),('side_minus_y',-y)]:
        first=None
        for distance in np.arange(.1,2.001,.1):
            moved=ball.translate(axis*distance);hits=[{'part':n,'overlap_mm3':float((moved^s).volume())} for n,s in fixed.items()];hits=[h for h in hits if h['overlap_mm3']>.001]
            if hits:first={'translation_mm':float(distance),'contacts':hits};break
        rows.append({'direction':name,'first_sampled_interference':first})
    report={'scope':'ball-only straight translations, fixed assembled socket/caps; not an all-path escape proof or retention-load test','step_mm':.1,'maximum_probe_mm':2,'directions':rows,'all_six_straight_directions_blocked_within_2mm':all(r['first_sampled_interference'] is not None for r in rows),'requires':'real cap fasteners, strength and printed clearance validation'};(DEST/'capture_probes.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)

if __name__=='__main__':run()
