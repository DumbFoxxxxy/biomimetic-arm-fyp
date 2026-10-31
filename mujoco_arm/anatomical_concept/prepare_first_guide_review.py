"""Check local fitting and both selected guide placements, export real scene."""
from screen_first_guides import *
from build_first_guide import write

def run():
    report=json.loads((DEST/'first_guide_placement_screen.json').read_text());parts={k:load(TEMPLATE/(k+'.obj'))[0] for k in json.loads((TEMPLATE/'design.json').read_text())['parts']}
    paths=[];fixed=['cassette_fork']
    for label,keys,direction in [('wheel_and_sleeve_side_entry',['grooved_wheel','shoulder_sleeve'],[0,1,0]),('rear_washer_and_loose_nut',['rear_washer','M3_nut'],[0,0,-1]),('front_washer_and_axle',['front_washer','M3x16_axle'],[0,0,1])]:
        first=None
        distances=np.r_[np.arange(20,5,-1),np.arange(5,-.00001,-.1)]
        for distance in distances:
            h=[]
            for key in keys:h.extend({'moving':key,**hit} for hit in hits(parts[key].translate(np.array(direction)*distance),{k:parts[k] for k in fixed}))
            if h:first={'approach_distance_mm':float(distance),'hits':h};break
        paths.append({'label':label,'moving_parts':keys,'direction':direction,'first_obstruction':first,'sample_count':len(distances) if first is None else None,'final_5mm_step_mm':.1})
        fixed+=keys
    (DEST/'cassette_assembly_checks.json').write_text(json.dumps({'scope':'ideal rigid local insertion paths; loose nut manually held, tools, cable installation, tolerances and frame absent','paths':paths,'all_clear':all(p['first_obstruction'] is None for p in paths)},indent=2),encoding='utf8')
    assert all(p['first_obstruction'] is None for p in paths)
    print('guide local fitting paths',len(paths),'clear',flush=True)
    selected={m['muscle']:m['best_placement_only'] for m in report['modules']};placed={};rows=[];geoms=[]
    folder=DEST/'paired_placement';folder.mkdir(exist_ok=True)
    for fid,r in selected.items():
        assert r and r['continuous_forward_prefix_deg']==150
        F=np.array(r['frame_columns']);Q=np.array(r['centre_display_mm']);T=np.c_[F,Q]
        for key,s in parts.items():placed[fid+'_'+key]=s.transform(T)
        wire,_,_,_=route(r['source_distance_mm'],r['offset_angle_deg'],r['outgoing_angle_deg']);placed[fid+'_wire_corridor']=wire.transform(T)
    cross=[]
    for (a,sa),(b,sb) in itertools.combinations(placed.items(),2):
        if a[:6]==b[:6]:continue
        v=intersection(sa,sb)
        if v>.001:cross.append({'pair':[a,b],'volume_mm3':v})
    for key,s in placed.items():
        rows.append(write(folder/(key+'.stl'),s))
        colour='.12 .15 .18 1' if key.endswith('grooved_wheel') else '.75 .77 .8 1'
        if key.endswith('wire_corridor'):colour='1 .25 .05 1'
        geoms.append((key,colour))
    assets=''.join(f'<mesh name="{k}" file="{k}.obj" scale=".001 .001 .001"/>' for k,_ in geoms)
    scene=''.join(f'<geom name="{k}" type="mesh" mesh="{k}" rgba="{c}" contype="0" conaffinity="0"/>' for k,c in geoms)
    (folder/'guides.xml').write_text(f'<mujoco model="M6 ungrounded guide placements"><compiler meshdir="."/><visual><global offwidth="1200" offheight="800"/></visual><asset>{assets}</asset><worldbody><light pos="0 -1 2" dir="0 1 -2" diffuse=".8 .8 .8"/>{scene}</worldbody></mujoco>',encoding='utf8')
    (folder/'selection.json').write_text(json.dumps({'scope':'two placed guides plus finite cable corridors only; no ground frame or downstream support, not manufactured assembly or actuation','placements':selected,'cross_module_interferences':cross,'frozen':False},indent=2),encoding='utf8')
    (folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[r for r in rows if not r['name'].endswith('wire_corridor')]}),encoding='utf8')
    print('selected guide cross-module hits',len(cross),flush=True);assert not cross

if __name__=='__main__':run()
