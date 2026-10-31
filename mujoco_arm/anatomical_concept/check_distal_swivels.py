"""M5 neutral, actual articulation, finite retention and assembly checks."""
from build_distal_swivels import *
from scipy.spatial.transform import Rotation
import sys

def run():
    assembly=json.loads((BASE/'assembly.json').read_text(encoding='utf8'));rows={r['key']:r for r in assembly['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};modules={}
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=json.loads((folder/'design.json').read_text(encoding='utf8'));parts={name:load(folder/(name+'.obj'))[0] for name in design['parts']};modules[fid]=(design,parts)
        del world[design['source_cap_key']]
        for key,s in parts.items():world[fid+'_swivel_'+key]=s
    hard=[];soft=[]
    for (a,sa),(b,sb) in itertools.combinations(world.items(),2):
        volume=intersection(sa,sb)
        if volume>.001:(soft if any(k in rows and rows[k]['kind']=='TPU' for k in [a,b]) else hard).append({'pair':[a,b],'volume_mm3':volume})
    neutral={'scope':'all131 manufactured neutral instances, including two floating endpoint assemblies; not installed whole-muscle motion','instance_count':len(world),'rigid_interferences':hard,'TPU_contact_review_flags':soft,'rigid_threshold_mm3':.001,'frozen':False}
    (DEST/'neutral_checks.json').write_text(json.dumps(neutral,indent=2),encoding='utf8');print('M5 neutral',len(world),'hard',len(hard),'soft',len(soft),flush=True)
    for fid,(design,parts) in modules.items():
        folder=DEST/fid;selection=design['selected'];C=np.array(selection['ball_centre_display_mm']);u=np.array(selection['neutral_terminal_axis']);v=np.array(selection['retainer_v_axis']);w=np.cross(u,v);terminal=parts['ball_crimp_blank'];others={k:s for k,s in world.items() if k!=fid+'_swivel_ball_crimp_blank'}
        cached=None
        if '--assembly-only' in sys.argv and (folder/'mechanism_checks.json').exists():
            cached=json.loads((folder/'mechanism_checks.json').read_text(encoding='utf8'))
            for key,digest in cached['source_mesh_sha256'].items():assert hashlib.sha256((folder/(key+'.obj')).read_bytes()).hexdigest()==digest
        def collision(s,against):
            hard=[];soft=[]
            for key,part in against.items():
                volume=intersection(s,part)
                if volume>.001:(soft if key in rows and rows[key]['kind']=='TPU' else hard).append({'part':key,'volume_mm3':volume})
            return hard,soft
        articulation=[]
        for azimuth in ([] if cached else range(0,360,15)):
            rotaxis=v*np.cos(np.radians(azimuth))+w*np.sin(np.radians(azimuth));path=[];first=None
            for angle in range(26):
                R=Rotation.from_rotvec(rotaxis*np.radians(angle)).as_matrix();s=terminal.transform(np.c_[R,C-R@C]);hits,contacts=collision(s,others);row={'angle_deg':angle,'rigid_hits':hits,'TPU_contact_flags':contacts};path.append(row)
                if hits and first is None:first=angle
            last=25. if first is None else 0.;boundary=None
            if first is not None and first>0:
                low,high=first-1.,float(first)
                while high-low>.1:
                    mid=(low+high)/2;R=Rotation.from_rotvec(rotaxis*np.radians(mid)).as_matrix();hits,contacts=collision(terminal.transform(np.c_[R,C-R@C]),others)
                    if hits:high=mid;boundary={'angle_deg':mid,'hits':hits}
                    else:low=mid
                last=low
            back=[]
            for angle in range(25,-1,-1):
                R=Rotation.from_rotvec(rotaxis*np.radians(angle)).as_matrix();hits,contacts=collision(terminal.transform(np.c_[R,C-R@C]),others);back.append({'angle_deg':angle,'rigid_hits':hits,'TPU_contact_flags':contacts})
            articulation.append({'azimuth_deg':azimuth,'continuous_clear_prefix_deg':last,'first_boundary':boundary,'forward':path,'return':back})
        # Finite straight pullout samples; not a proof for all escape curves.
        directions=[unit(np.array(p,dtype=float)) for p in itertools.product([-1,0,1],repeat=3) if any(p)]+[u,-u]
        retention=[];housing={name:s for name,s in parts.items() if name!='ball_crimp_blank'}
        for direction in ([] if cached else directions):
            first=None
            for offset in np.arange(.1,3.001,.1):
                hits,contacts=collision(terminal.translate(direction*offset),housing)
                if hits:first={'offset_mm':float(offset),'hits':hits};break
            retention.append({'direction':direction.tolist(),'first_obstruction':first})
        if cached:articulation=cached['articulation_paths'];retention=cached['retention_directions']
        clamp_stops=[]
        if design.get('integral_clamp_stops',{}).get('present'):
            for i,point in enumerate(selection['M2_centres_display_mm']):
                p=np.array(point);stop=rod(p-u*.15,p+u*.15,2.5)-rod(p-u,p+u,1.1);missing=max(0,float((stop-parts['distal_socket_base']).volume()));gap=float(stop.min_gap(parts['socket_retainer'],1));close=intersection(stop,parts['socket_retainer'].translate(-u*.01))
                clamp_stops.append({'ear':i,'missing_stop_material_mm3':missing,'nominal_stop_to_retainer_gap_mm':gap,'lid_0_01mm_closing_overlap_mm3':close,'pass':missing<.001 and gap<1e-4 and close>.1})
        module={key:s for key,s in world.items() if key.startswith(fid+'_') and '_bone_' not in key and '_swivel_' not in key};all_parts={**module,**parts};soft_names={key for key in module if key in rows and rows[key]['kind']=='TPU'}
        def probe(label,moving,fixed,direction,distance):
            samples=[]
            offsets=np.r_[np.arange(distance,5,-1),np.linspace(min(5,distance),0,int(min(5,distance)*10)+1)]
            for offset in offsets:
                hits=[]
                for a in moving:
                    moved=all_parts[a].translate(direction*offset)
                    for b in fixed:
                        volume=intersection(moved,all_parts[b])
                        if volume>.001:hits.append({'pair':[a,b],'volume_mm3':volume,'classification':'soft_review' if b in soft_names else 'rigid'})
                samples.append({'offset_mm':float(offset),'hits':hits})
            return {'label':label,'moving':moving,'direction':direction.tolist(),'samples':samples,'rigid_clear':not any(h['classification']=='rigid' for r in samples for h in r['hits']),'soft_flags':sum(h['classification']=='soft_review' for r in samples for h in r['hits']),'first_rigid_limit':next((r for r in samples if any(h['classification']=='rigid' for h in r['hits'])),None)}
        source=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));z=np.array(source['basis_columns_display'])[:,2]
        fixed=[k for k in module if '_distal_bolt_' not in k and '_distal_head_washer_' not in k and not k.endswith('active_wire_crimp_barrel')];paths=[probe('modified_cap_to_TPU_before_head_fasteners',['distal_socket_base'],fixed,z,25)]
        fixed+=['distal_socket_base']
        for side in ['left','right']:
            pair=[fid+'_distal_bolt_'+side,fid+'_distal_head_washer_'+side];paths.append(probe('original_M3_'+side,pair,fixed,z,25));fixed+=pair
        key=fid+'_active_wire_crimp_barrel';paths.append(probe('original_active_wire_barrel',[key],fixed,z,10));fixed.append(key)
        paths.append(probe('ball_before_retainer',['ball_crimp_blank'],fixed,u,20));fixed.append('ball_crimp_blank')
        paths.append(probe('socket_retainer',['socket_retainer'],fixed,u,20));fixed.append('socket_retainer')
        for side in ['left','right']:
            # Position loose rear hardware before a shaft would block lateral
            # access. These pieces need manual holding; no captive nut implied.
            pair=[side+'_nut',side+'_back_washer'];options=[]
            sign=-1 if side=='left' else 1
            for label,direction in [('axial',-u),('lateral_v',v*sign),('lateral_w',w),('minus_w',-w),('minus_v',-v*sign),('diagonal',unit(v*sign+w)),('other_diagonal',unit(v*sign-w))]:
                option=probe('M2_nut_'+side+'_'+label,pair,fixed,direction,20);options.append(option)
                if option['rigid_clear'] and not option['soft_flags']:break
            chosen=min(options,key=lambda r:(not r['rigid_clear'],r['soft_flags']))
            chosen['tested_access_directions']=[{'label':r['label'],'rigid_clear':r['rigid_clear'],'soft_flags':r['soft_flags'],'first_rigid_limit':r['first_rigid_limit']} for r in options]
            paths.append(chosen);fixed+=pair
            pair=[side+'_M2x10',side+'_head_washer'];paths.append(probe('M2_bolt_'+side,pair,fixed,u,20));fixed+=pair
        report={'scope':'selected real endpoint:neutral assembly,24 angular radial paths and28 finite straight pullout directions; not strength, escape-proofing, threading, crimp or compliant arm actuation','articulation_paths':articulation,'angle_step_deg':1,'boundary_threshold_deg':.1,'minimum_sampled_articulation_prefix_deg':min(p['continuous_clear_prefix_deg'] for p in articulation),'retention_directions':retention,'retention_step_mm':.1,'all_tested_pullout_directions_obstructed_by_3mm':all(r['first_obstruction'] for r in retention),'assembly_paths':paths,'assembly_final5mm_step_mm':.1,'assembly_other_step_mm':1,'all_rigid_assembly_paths_clear':all(p['rigid_clear'] for p in paths),'assembly_soft_flags':sum(p['soft_flags'] for p in paths),'source_mesh_sha256':{key:hashlib.sha256((folder/(key+'.obj')).read_bytes()).hexdigest() for key in parts},'frozen':False}
        report['integral_clamp_stop_checks']=clamp_stops
        report['all_integral_clamp_stops_pass']=bool(clamp_stops) and all(p['pass'] for p in clamp_stops)
        (folder/'mechanism_checks.json').write_text(json.dumps(report,indent=2),encoding='utf8')
        print(fid,'articulation',report['minimum_sampled_articulation_prefix_deg'],'retention',report['all_tested_pullout_directions_obstructed_by_3mm'],'fitting',sum(p['rigid_clear'] for p in paths),'/',len(paths),'soft',report['assembly_soft_flags'],flush=True)
        for p in paths:
            if not p['rigid_clear']:print(p['label'],p['first_rigid_limit'],flush=True)
if __name__=='__main__':run()
