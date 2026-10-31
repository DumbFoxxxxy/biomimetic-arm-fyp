"""M7 integral humeral-adapter extensions and keyed first-sheave receivers.

Real material grounds each fork to an existing two-bolt humeral adapter; no
new bone holes or joint-centre moves. Strength, cap motion and drive unproven.
"""
from screen_first_guides import *
from build_first_guide import annular_profile, write
import sys
SWEPT_ROOT='--swept-root' in sys.argv
LOW_ROOT='--low-root' in sys.argv or SWEPT_ROOT
DEST7=DEST/('ground_frame_swept' if SWEPT_ROOT else 'ground_frame_low_root' if LOW_ROOT else 'ground_frame_candidate')

def stepped_mount_bolt(length=20):
    n=32;vertices=[];faces=[]
    for r,z in [(1.5,5.5-length),(1.5,5.5),(2.75,5.5),(2.75,8.5)]:
        vertices.extend([r*np.cos(a),r*np.sin(a)-13,z] for a in np.arange(n)*2*np.pi/n)
    vertices.extend([[0,-13,5.5-length],[0,-13,8.5]])
    for i in range(3):
        for j in range(n):
            a=i*n+j;b=i*n+(j+1)%n;c=b+n;d=a+n;faces.extend([[a,b,c],[a,c,d]])
    for j in range(n):faces.extend([[(j+1)%n,j,4*n],[3*n+j,3*n+(j+1)%n,4*n+1]])
    s=mf.Manifold(mf.Mesh64(np.array(vertices),np.array(faces,dtype=np.uint64)));assert s.status()==mf.Error.NoError;return s

def local_receiver():
    base=mf.Manifold.cube([11,8,3]).translate([-5.5,-17.5,-8])
    left=mf.Manifold.cube([2.35,8,3]).translate([-5.5,-17.5,-5])
    right=mf.Manifold.cube([2.35,8,3]).translate([3.15,-17.5,-5])
    rear=mf.Manifold.cube([11,1.35,3]).translate([-5.5,-17.5,-5])
    # Open +Y insertion path; rear plane contacts key at Z=-5 exactly.
    receiver=(base+left+right+rear-rod([0,-13,-12],[0,-13,0],1.7,32)).as_original().simplify(.001)
    hardware={'key_M3x20':stepped_mount_bolt(),
              'key_front_washer':annular_profile([(3.5,5),(3.5,5.5),(1.6,5.5),(1.6,5)],32).translate([0,-13,0]),
              'key_back_washer':annular_profile([(3.5,-8.5),(3.5,-8),(1.6,-8),(1.6,-8.5)],32).translate([0,-13,0]),
              'key_nut':mf.Manifold.cylinder(2.4,5.5/np.sqrt(3),5.5/np.sqrt(3),6).translate([0,-13,-10.9])-rod([0,-13,-12],[0,-13,-7],1.6,32)}
    if LOW_ROOT:
        # Rear-entering shaft; front nut/washer can be placed laterally first.
        R=np.diag([1.,-1.,-1.]);T=np.c_[R,[0,-26,-3]]
        hardware.pop('key_M3x20')
        hardware['key_M3x18']=stepped_mount_bolt(18).transform(T)
        hardware['key_nut']=hardware['key_nut'].transform(T)
    return receiver,hardware

def path_solid(points,r=3):
    return mf.Manifold.batch_boolean([rod(a,b,r) for a,b in zip(points,points[1:])]+[mf.Manifold.sphere(r,16).translate(p) for p in points[1:-1]],mf.OpType.Add)

def run():
    DEST7.mkdir(exist_ok=True);assembly=json.loads((CURRENT/'assembly.json').read_text());selection=json.loads((DEST/'paired_placement/selection.json').read_text())['placements']
    world={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in assembly['instances']};kinds={r['key']:r['kind'] for r in assembly['instances']};placements={r['id']:r for r in assembly['placements']}
    # Include both placed cassettes and actual finite neutral cable corridors.
    for p in (DEST/'paired_placement').glob('*.obj'):world[p.stem]=load(p)[0];kinds[p.stem]='cable_corridor' if p.stem.endswith('wire_corridor') else 'rigid'
    support=json.loads((ROOT/'骨骼/支撑参考/sources.json').read_text(encoding='utf8'))
    for r in support['models']:
        if r['kind']=='bone_reference' and r['id'] not in world:
            world['reference_'+r['id']]=load(OUT/'repaired_supports'/(r['id']+'.obj'),raw=True)[0];kinds['reference_'+r['id']]='bone_reference'
    receiver,hardware=local_receiver();reports=[];chosen={};rows=[]
    definition=json.loads((OUT/'joint_model/definition.json').read_text());shift=np.array(definition['source_coordinate_shift_mm']);G=np.array(definition['centres_source_mm']['GH'])-shift;E=np.array(definition['centres_source_mm']['elbow'])-shift;long=unit(E-G)
    axis=np.array(json.loads((OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection/assembly.json').read_text())['axis_source'])
    moving={r['key']:world[r['key']] for r in assembly['instances'] if r['owner']=='elbow' and not r['key'].endswith('_tendon_stop')}
    static_shoulder={k:s for k,s in world.items() if k.startswith('reference_') or any(r['key']==k and r['owner'] in ['world','SC','AC'] for r in assembly['instances'])}
    def axis_at_z(z):return G+long*((z-G[2])/long[2])
    for fid in ['FJ1486','FJ1477']:
        r=selection[fid];F=np.array(r['frame_columns']);Q=np.array(r['centre_display_mm']);T=np.c_[F,Q];oldkey=fid+'_bone_adapter';original=world[oldkey]
        # Existing outboard brace is mechanically connected to BOTH bone ears.
        source_path=np.array(placements[fid]['brace_paths_display_mm'][1])
        roots=np.array([source_path[-2],source_path[-2]*.75+source_path[-1]*.25]) if LOW_ROOT else source_path[1:3]
        side=1 if fid=='FJ1486' else -1
        socket=receiver.transform(T);hardware_world={fid+'_'+k:s.transform(T) for k,s in hardware.items()}
        for key,s in hardware_world.items():
            for hit in hits(s,{k:b for k,b in world.items() if k!=oldkey}):
                raise RuntimeError('Key mount hardware obstructed '+key+' '+str(hit))
        trials=[];best=None
        parameters=itertools.product([20,30,40,50],[25,35,45],[0,20,-20]) if LOW_ROOT else itertools.product([0,20,40,60],[0,20,40],[0,20,-20])
        for outx,offsety,lower_offset in parameters:
            # Receiver attaches outboard of the key, not through its bolt/nut.
            end=Q+F@np.array([side*9,-13,-6.5]);socket_root=Q+F@np.array([side*4.7,-13,-6.5])
            # A front-side longitudinal rail crosses the lever's swept annulus.
            # Route the candidate through posterior free space instead; actual
            # solids must pass the elbow and GH gates before being accepted.
            y=(1 if SWEPT_ROOT or fid=='FJ1477' else -1)*(offsety if LOW_ROOT else 60+offsety)
            upper_z=float(roots[:,2].max()+(0 if LOW_ROOT else 10));lower_z=float(end[2]+lower_offset)
            upper=axis_at_z(upper_z)+[outx,y,0] if LOW_ROOT else np.array([outx,y,upper_z])
            lower=axis_at_z(lower_z)+[outx,y,0] if LOW_ROOT else np.array([outx,y,lower_z])
            paths=[np.vstack([root,upper,lower,end]) for root in roots]
            extension=mf.Manifold.batch_boolean([socket,rod(socket_root,end,2.5)]+[path_solid(p) for p in paths],mf.OpType.Add)
            combined=(original+extension).as_original().simplify(.001)
            components=len([s for s in combined.decompose() if s.volume()>.001]);collisions=hits(extension,{k:b for k,b in world.items() if k!=oldkey})
            collisions+=hits(extension,hardware_world)
            row={'outboard_x_mm':outx,'y_offset_parameter_mm':offsety,'lower_z_offset_mm':lower_offset,'components':components,'interferences':collisions,'paths_display_mm':[p.tolist() for p in paths],'receiver_root_display_mm':socket_root.tolist(),'extension_volume_mm3':extension.volume()}
            if SWEPT_ROOT and components==1 and not collisions:
                fixed_new=mf.Manifold.batch_boolean([extension]+list(hardware_world.values()),mf.OpType.Add);limits={}
                for angle in range(151):
                    R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();h=hits(fixed_new,{k:s.transform(np.c_[R,E-R@E]) for k,s in moving.items()})
                    if h:limits['elbow']={'angle_deg':angle,'hits':h};break
                if not limits:
                    for label,angles in [('GH_positive',list(range(83))+[82.25]),('GH_negative',range(0,-91,-1))]:
                        for angle in angles:
                            R=Rotation.from_rotvec(long*np.radians(angle)).as_matrix();h=hits(fixed_new.transform(np.c_[R,G-R@G]),static_shoulder)
                            if h:limits[label]={'angle_deg':angle,'hits':h};break
                row['new_solids_motion_screen_limits']=limits
                (DEST7/(fid+'_search_checkpoint.json')).write_text(json.dumps({'trials':trials+[row],'scope':'new solids elbow and GH axial gate only, before complete suite'},indent=2),encoding='utf8')
                if limits:trials.append(row);print(fid,'trial',len(trials),'motion rejected',list(limits),flush=True);continue
            trials.append(row)
            if components==1 and not collisions:
                best=(row,combined,extension,hardware_world);break
        reports.append({'muscle':fid,'trials':trials,'selected':best[0] if best else None})
        print(fid,'frame trials',len(trials),'clear',bool(best),flush=True)
        if best:
            row,combined,extension,hardware_world=best;chosen[oldkey]=combined
            rows.append(write(DEST7/(oldkey+'.stl'),combined));rows.append(write(DEST7/(fid+'_extension_only.stl'),extension))
            for key,s in hardware_world.items():chosen[key]=s;rows.append(write(DEST7/(key+'.stl'),s))
            # Subsequent module must clear this actual chosen frame/hardware.
            world[oldkey]=combined;world.update(hardware_world)
    report={'scope':'M7 integral adapter extensions to actual keyed guide receivers; geometry candidate only, no active muscle motion/strength claim','modules':reports,'source_assembly_sha256':digest(CURRENT/'assembly.json'),'joint_centres_changed':False,'new_bone_holes':False,'guide_key_back_plane_local_z_mm':-5,'keyed_side_clearance_each_mm':.15,'nominal_rail_radius_mm':3,'mount_bolt':'M3x20 nominal envelope, real thread/drive not modeled','full_transmission_exists':False,'frozen':False}
    report['root_revision']='lower brace only, rear-entering bolt/front loose nut' if LOW_ROOT else 'upper/lower brace, front-entering bolt/rear nut'
    if LOW_ROOT:report['mount_bolt']='M3 nominal18mm shaft, rear entry, front loose nut; supplier/thread/drive unselected'
    (DEST7/'design.json').write_text(json.dumps(report,indent=2),encoding='utf8');(DEST7/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[r for r in rows if not r['name'].endswith('_extension_only')]}),encoding='utf8')
    if len(chosen)==10:
        mutual=[]
        for (ka,sa),(kb,sb) in itertools.combinations(chosen.items(),2):
            v=intersection(sa,sb)
            if v>.001:mutual.append({'pair':[ka,kb],'volume_mm3':v})
        (DEST7/'mutual_checks.json').write_text(json.dumps({'new_piece_count':10,'interferences':mutual},indent=2),encoding='utf8');print('M7 mutual hits',len(mutual),flush=True)

if __name__=='__main__':run()
