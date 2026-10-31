"""Finite rigid fitting paths for the actual end components, not threaded mechanics."""
from parameterize_elbow_muscles import *

def run(muscles=('FJ1486','FJ1477')):
    for fid in muscles:
        folder=DEST/fid;data=json.loads((folder/'design.json').read_text(encoding='utf8'));F=np.array(data['basis_columns_display']);parts={}
        for name in data['components']:
            v,f=read_obj(str(folder/(name+'.obj')));parts[name]=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)))
        def probe(moving,fixed,direction,distance):
            samples=[]
            for offset in range(distance,-1,-1):
                hits=[]
                for a in moving:
                    s=parts[a].translate(direction*offset)
                    for b in fixed:
                        overlap=max(0,float((s^parts[b]).volume()))
                        if overlap>.001:hits.append({'parts':[a,b],'volume_mm3':overlap})
                samples.append({'offset_mm':offset,'hits':hits})
            return {'moving':moving,'direction_display':direction.tolist(),'step_mm':1,'samples':samples,'collision_free':not any(s['hits'] for s in samples)}
        rows=[]
        for end,sign in [('proximal',-1),('distal',1)]:
            cap=end+'_rigid_cap'
            fixed=[n for n in parts if n!=cap and '_bolt_' not in n and '_head_washer_' not in n and n!='active_wire_crimp_barrel']
            rows.append(probe([cap],fixed,F[:,2]*sign,20))
        for end,sign in [('proximal',-1),('distal',1)]:
            for side,xsign in [('left',-1),('right',1)]:
                moving=[end+'_nut_'+side,end+'_nut_washer_'+side]
                fixed=[n for n in parts if n not in moving and '_bolt_' not in n]
                rows.append(probe(moving,fixed,F[:,0]*xsign,20))
                bolt=end+'_bolt_'+side;rows.append(probe([bolt],[n for n in parts if n!=bolt],F[:,2]*sign,20))
        rows.append(probe(['active_wire_crimp_barrel'],[n for n in parts if n!='active_wire_crimp_barrel'],F[:,2],10))
        result={'scope':'1mm sampled ideal cap/nut/bolt/barrel fitting; thread engagement, crimp deformation, tool access and preload not simulated','paths':rows,'all_paths_collision_free':all(r['collision_free'] for r in rows),'assembly_order':['slide M3 nuts+inner washers through TPU access slots','seat proximal and distal caps before head washers/bolts/barrel','place head washers and insert bolts; nominal threads only','thread active wire, crimp barrel with a real tool and seat it','fit passive tendon through distal transverse eyelet; external path pending'],'frozen':False}
        (folder/'end_assembly_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(fid,'fitting paths',sum(r['collision_free'] for r in rows),'/',len(rows),'failed',[r['moving'] for r in rows if not r['collision_free']],flush=True)

if __name__=='__main__':run()
