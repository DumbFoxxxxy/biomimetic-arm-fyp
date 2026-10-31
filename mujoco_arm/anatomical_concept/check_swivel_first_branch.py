"""Check a straight30mm output corridor; does not invent a turning guide."""
from build_distal_swivels import *
def run():
    assembly=json.loads((DEST/'assembly.json').read_text(encoding='utf8'));parts={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in assembly['instances']};kinds={r['key']:r['kind'] for r in assembly['instances']};reports=[]
    for fid in ['FJ1486','FJ1477']:
        design=json.loads((DEST/fid/'design.json').read_text(encoding='utf8'));C=np.array(design['selected']['ball_centre_display_mm']);u=np.array(design['selected']['neutral_terminal_axis']);active=np.array(design['neutral_active_wire_pull_direction']);wire=rod(C+u*3.5,C+u*30,.3);hard=[];soft=[]
        for key,s in parts.items():
            volume=intersection(wire,s)
            if volume>.001:(soft if kinds[key]=='TPU' else hard).append({'part':key,'volume_mm3':volume})
        reports.append({'muscle':fid,'neutral_output_axis':u.tolist(),'opposition_to_neutral_active_pull':-float(u@active),'corridor_endpoint_display_mm':(C+u*30).tolist(),'checked_wire_span_from_ball_centre_mm':[3.5,30],'rigid_interferences':hard,'static_TPU_contact_flags':soft,'actual_turning_guide_exists':False,'external_force_or_grip_verified':False})
    (DEST/'first_branch_checks.json').write_text(json.dumps({'scope':'neutral straight output corridor and traction-direction projection only; endpoint is not a support/control point, no guide constraint or TPU law imposed','modules':reports,'assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest(),'frozen':False},indent=2),encoding='utf8');print(reports)
if __name__=='__main__':run()
