"""M7 neutral tangent-aligned free terminals and first supported cable branches.

CAD neutral packing retains nominal M5B ball orientation; this records the
explicit movable-ball pose used by the static routing illustration.
"""
from build_guide_ground_frames import *

def run():
    a=json.loads((DEST7/'assembly.json').read_text());world={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in a['instances']};kinds={r['key']:r['kind'] for r in a['instances']};selection=json.loads((DEST/'paired_placement/selection.json').read_text())['placements'];poses=[]
    for fid,s in selection.items():
        d=json.loads((CURRENT/fid/'design.json').read_text());C=np.array(d['selected']['ball_centre_display_mm']);u=np.array(d['selected']['neutral_terminal_axis']);F=np.array(s['frame_columns']);_,inc,_,_=route(s['source_distance_mm'],s['offset_angle_deg'],s['outgoing_angle_deg']);direction=F@inc;v=np.cross(u,direction);K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]);R=np.eye(3)+K+K@K/(1+u@direction);key=fid+'_swivel_ball_crimp_blank';world[key]=world[key].transform(np.c_[R,C-R@C]);poses.append({'muscle':fid,'part':key,'rotation_matrix':R.tolist(),'pivot_display_mm':C.tolist(),'nominal_to_tangent_angle_deg':s['terminal_tilt_deg']})
    checks=[]
    for pose in poses:
        key=pose['part'];h=hits(world[key],{k:s for k,s in world.items() if k!=key});checks.append({'part':key,'interferences':h})
    for fid in selection:
        p=DEST/'paired_placement'/(fid+'_wire_corridor.obj');s=load(p)[0];checks.append({'part':fid+'_neutral_supported_branch','interferences':hits(s,world)})
    result={'scope':'two static neutral free-ball tangent poses and first cable branch against all153 manufactured parts, actual grounded guide bends; outgoing endpoints free, no full cap/muscle motion','checks':checks,'poses':poses,'all_clear':not any(r['interferences'] for r in checks),'CAD_neutral_packing_ball_pose':'original M5B neutral, same geometry; renderer explicitly articulates free balls to these tangent poses','source_assembly_sha256':digest(DEST7/'assembly.json'),'source_sha256':{str(ROOT/r['file']):digest(ROOT/r['file']) for r in a['instances']},'frozen':False};(DEST7/'neutral_supported_branch_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M7 neutral tangent poses/branches',len(checks),'clear',result['all_clear'],flush=True);assert result['all_clear']

if __name__=='__main__':run()
