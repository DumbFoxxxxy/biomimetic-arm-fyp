"""Search a real curved cap outlet for the FJ1467 external lead.

Preserves the bone seat and anatomy. Circular bend radii are design candidates,
not certified steel or conduit limits. Does not deploy a manufacturing baseline.
"""
import copy
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import build_shoulder_tendon_segments as T
import shoulder_pads
import interference as I
import check_shoulder_assembly as C
from build_shoulder_unified_muscle import local_export


def curved_outlet(inner, normal, lateral, thickness, tilt_deg, bend_radius):
    """Tangent-continuous circular bend from the existing axial inlet.

    Curve is sampled at <=1 degree and explicitly crosses both cap faces.
    """
    theta=np.radians(tilt_deg)
    t=bend_radius*np.tan(theta/2)
    # There must be straight material below the existing 2mm inlet stub and
    # a full bend before the outlet crosses the top cap face.
    if t>=2 or t*np.cos(theta)>thickness:return None
    direction=normal*np.cos(theta)+lateral*np.sin(theta)
    alpha_cross=np.arcsin(t/bend_radius)
    alpha=np.unique(np.r_[np.linspace(0,theta,int(np.ceil(tilt_deg))+1),alpha_cross])
    centre=inner-normal*t+lateral*bend_radius
    arc=centre+np.sin(alpha)[:,None]*normal*bend_radius-np.cos(alpha)[:,None]*lateral*bend_radius
    index=int(np.argmin(abs(alpha-alpha_cross)))
    outer=inner+direction*thickness/np.cos(theta)
    end=outer+direction*5
    guide=np.vstack([arc,outer,end])
    return dict(guide=guide,arc=arc,cross_index=index,outer=outer,end=end,
                inner=arc[index],direction=direction)


def main():
    source=D.DEST/'routes_separated.json'
    spec=json.loads(source.read_text(encoding='utf8'))
    row=next(r for r in spec['routes'] if r['head']=='FJ1467')
    old=row['ports']['origin'];base=np.array(old['inner_face_source_mm'])
    n=np.array(old['normal_source']);basis=shoulder_pads.pad_basis(n)
    thickness=np.dot(np.array(old['outer_face_source_mm'])-base,n)
    bones={k:D.SP.load_solid(T.BONE_DIR/v) for k,v in D.SP.BONES.items()}
    trials=[];selected=None
    # R3 is a declared geometric candidate only; larger bends may require a
    # different inlet stub. Do not silently certify the unknown cable material.
    for tilt in [30,45,60]:
        for phase in range(0,360,15):
            a=np.radians(phase);u=np.cos(a)*basis[:,0]+np.sin(a)*basis[:,1]
            choice=curved_outlet(base,n,u,thickness,tilt,3.)
            if choice is None:continue
            guide=D.outer_envelope(choice['guide'],1.8)
            hits={k:I.intersection_volume(guide,b) for k,b in bones.items()}
            trial=dict(tilt_deg=tilt,phase_deg=phase,bend_radius_mm=3.,
                       outer_radius_mm=1.8,bone_intersections_mm3=hits)
            trials.append(trial)
            if all(v<=.001 for v in hits.values()):
                if selected is None or sum(hits.values())<sum(selected[2]['bone_intersections_mm3'].values()):
                    selected=(choice,guide,trial)
                if not any(hits.values()):break
        if selected and not any(selected[2]['bone_intersections_mm3'].values()):break
    folder=D.DEST/'origin_outlet_candidate';folder.mkdir(exist_ok=True)
    report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                scope='curved rigid-cap outlet candidate; not complete TPU or installed conduit',
                trials=trials,selected=None,stage_pass=False)
    if selected:
        choice,guide,trial=selected
        updated=copy.deepcopy(spec)
        new=next(r for r in updated['routes'] if r['head']=='FJ1467')
        # Split the circular bend at the actual inner cap face. Retain the
        # existing inward 2mm stub and every subsequent muscle/core point.
        arc=choice['arc'];j=choice['cross_index']
        connector=np.vstack([arc[:j+1][::-1],np.array(row['connectors']['origin']['path_source_mm'])[1:]])
        # Wire continues 3mm past the guide body, so the lumen is genuinely
        # open at its rounded end instead of making another sealed tube.
        free_end=choice['end']+choice['direction']*3
        full=np.vstack([free_end,choice['end'],choice['outer'],arc[j+1:][::-1],connector,
                        np.array(row['core_path_source_mm'])[1:],
                        np.array(row['connectors']['distal']['path_source_mm'])[-2::-1],
                        row['ports']['distal']['outer_face_source_mm']])
        new['ports']['origin'].update(inner_face_source_mm=choice['inner'].tolist(),
            offset_from_seat_mm=(choice['inner']-base).tolist(),
            outer_face_source_mm=choice['outer'].tolist(),hand_pull_lead_source_mm=free_end.tolist())
        new['connectors']['origin']['path_source_mm']=connector.tolist()
        new['connectors']['origin']['outer_envelope_bone_intersection_mm3']={
            k:I.intersection_volume(D.outer_envelope(connector,3.2),b) for k,b in bones.items()}
        new['connectors']['origin']['bone_clear_at_nominal_geometry']=all(
            v<=.001 for v in new['connectors']['origin']['outer_envelope_bone_intersection_mm3'].values())
        new['full_wire_path_source_mm']=full.tolist()
        new['rigid_origin_outlet']=dict(path_source_mm=choice['guide'].tolist(),
            outer_radius_mm=1.8,lumen_radius_mm=.86,bend_radius_mm=3.,material_verified=False)
        new['status']='candidate_requires_full_TPU_regeneration_and_external_reaction_mount'
        updated['scope']='static candidate with curved FJ1467 rigid origin outlet; not stage acceptance'
        updated['parent_routes_sha256']=report['source_sha256']
        deltoid=[r for r in updated['routes'] if r['module']=='deltoid']
        updated['deltoid_wire_pair_minimum_mm']=[dict(a=a['head'],b=b['head'],
            distance_mm=D.path_distance(np.array(a['full_wire_path_source_mm']),np.array(b['full_wire_path_source_mm'])))
            for a,b in itertools.combinations(deltoid,2)]
        (D.DEST/'routes_outlet_candidate.json').write_text(json.dumps(updated,indent=2),encoding='utf8')
        cap=C.load_stl(D.frames.OUTPUT_DIR/'shoulder_parts/FJ1467_origin_cap.stl')
        lumen=D.outer_envelope(full,.86)
        # The guide is actual material, continuous with the rigid cap, with
        # a matching drilled lumen. It cannot be assumed clear of unregenerated TPU.
        cap=(cap+guide)-lumen
        before_bone_trim={k:I.intersection_volume(cap,b) for k,b in bones.items()}
        for b in bones.values():cap=cap-b
        restored,meta=local_export(cap,folder/'FJ1467_origin_cap_with_guide.stl',base)
        steel=D.outer_envelope(full,.3)
        report['selected']=dict(**trial,part=meta,
            steel_cap_intersection_mm3=I.intersection_volume(steel,restored),
            cap_bone_trim_volume_mm3=before_bone_trim,
            steel_bone_intersections_mm3={k:I.intersection_volume(steel,b) for k,b in bones.items()},
            cap_bone_intersections_mm3={k:I.intersection_volume(restored,b) for k,b in bones.items()},
            guide_path_source_mm=choice['guide'].tolist())
    report['missing']=['TPU socket regenerated around actual guide outer body',
        'all other muscles, fasteners and installation clearance',
        'full route pair separation after modification',
        'steel/conduit bend qualification and friction',
        'conduit reaction mount, physical deformation and motion']
    (folder/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    selected_report=report['selected']
    print('Curved outlet selected:',{k:v for k,v in selected_report.items() if k!='guide_path_source_mm'} if selected_report else None,flush=True)


if __name__=='__main__':main()
