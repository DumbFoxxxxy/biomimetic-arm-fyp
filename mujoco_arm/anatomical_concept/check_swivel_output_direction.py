"""Work-sign diagnostic from actual neutral cable directions, not a rail law."""
from build_distal_swivels import *
def run():
    path=DEST/'straight_route_screen.json';data=json.loads(path.read_text(encoding='utf8'));assert data['assembly_sha256']==hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest()
    for summary in data['muscles']:
        fid=summary['muscle'];module=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));nodes=np.array(module['wire_polyline_display_mm']);active_pull=unit(nodes[-2]-nodes[-1]);projections=[]
        for row in data['samples']:
            case=next(r for r in row['muscles'] if r['muscle']==fid);direction=unit(np.array(case['actual_stop_reference_display_mm'])-np.array(case['ball_centre_display_mm']));opposition=-float(direction@active_pull)
            case['external_tension_opposition_to_neutral_active_pull']=opposition;projections.append(opposition)
        opposing=np.array(projections)>1e-6
        summary['neutral_active_pull_direction']=active_pull.tolist();summary['output_work_sign_diagnostic']={'definition':'opposition=-unit(B-C) dot neutral_active_terminal_pull; positive opposes nominal shortening, negative assists it. No force magnitude, imposed cap velocity or TPU motion proof.','opposition_at_elbow_zero':projections[0],'minimum':float(min(projections)),'maximum':float(max(projections)),'opposing_samples':int(opposing.sum()),'assisting_or_transverse_samples':int((~opposing).sum()),'first_nonopposing_elbow_deg':int(np.flatnonzero(~opposing)[0]) if not opposing.all() else None}
        print(fid,summary['output_work_sign_diagnostic'],flush=True)
    data['work_sign_scope']='Projection on the source-neutral last active cable span only. Does not establish deformed muscle force, torque capability or a motion constraint.'
    path.write_text(json.dumps(data,indent=2),encoding='utf8')
if __name__=='__main__':run()
