"""Revise the actual D17 fixture against the four integrated TPU modules.

Preserves posts and all bone-side positions. Joins intended structural unions
into two actual manufacturing solids, and relocates only the lateral base.
No claim about clamp force, stiffness, installation or physical retention.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import check_shoulder_assembly as C
import interference as I
from build_shoulder_unified_muscle import local_export


def main():
    parent=D.DEST/'whole_shoulder_static/actual_assembly_check.json'
    manifest=json.loads(parent.read_text(encoding='utf8'))
    solids={}
    for key,row in manifest['sources'].items():
        source=D.frames.PROJECT_ROOT/row['file']
        if hashlib.sha256(source.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Changed assembly input '+key)
        s=D.SP.load_solid(source) if source.suffix=='.obj' else C.load_stl(source)
        if row['translation_source_mm']:s=s.translate(row['translation_source_mm'])
        solids[key]=s
    fixtures={k.split(':',1)[1]:v for k,v in solids.items() if k.startswith('fixture:')}
    obstacles={k:v for k,v in solids.items() if not k.startswith('fixture:')}
    medial=fixtures['medial_base']+fixtures['medial_post']+fixtures['medial_beam']
    trials=[];chosen=None
    # Original export anchors cube's minimum corner at the post foot (not its
    # centre). Moving along -X keeps that foot inside the same full-size base.
    for displacement in [np.array([-float(x),0,0]) for x in range(21)]:
        base=fixtures['lateral_base'].translate(displacement)
        group=base+fixtures['lateral_post']
        connected=sum(p.volume()>1e-9 for p in group.decompose())==1
        hits={k:I.intersection_volume(base,v) for k,v in obstacles.items()}
        record=dict(displacement_source_mm=displacement.tolist(),connected=connected,
                    base_intersections_mm3=hits,base_minimum_gap_capped_0_5mm=None)
        trials.append(record)
        if not connected or any(v>.001 for v in hits.values()):continue
        gaps={k:I.gap(base,v,.5) for k,v in obstacles.items()}
        record['base_minimum_gap_capped_0_5mm']=gaps
        if min(gaps.values())<.5-1e-9:continue
        chosen=(group,record);break
    folder=D.DEST/'fixture_revision';folder.mkdir(exist_ok=True)
    report=dict(scope='two united rigid fixture groups with lateral base clearance; not clamp or strength verification',
        parent_report=str(parent.relative_to(D.frames.PROJECT_ROOT)),
        parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),trials=trials,selected=None,parts={},stage_pass=False)
    if chosen:
        lateral,selected=chosen;report['selected']=selected
        groups={'medial_frame':medial,'lateral_frame':lateral}
        restored={}
        for name,solid in groups.items():
            bb=np.array(solid.bounding_box());origin=(bb[:3]+bb[3:])/2
            restored[name],report['parts'][name]=local_export(solid,folder/(name+'.stl'),origin)
        report['actual_frame_material_intersections_mm3']={name:{k:I.intersection_volume(s,v) for k,v in obstacles.items()}
            for name,s in restored.items()}
        report['frame_pair_intersection_mm3']=I.intersection_volume(*restored.values())
        report['source_united_connections']=[dict(a=a,b=b,source_intersection_mm3=I.intersection_volume(fixtures[a],fixtures[b]))
            for a,b in [('lateral_base','lateral_post'),('medial_base','medial_post'),('medial_beam','medial_post')]]
        report['static_nonfixture_clear_at_0_001mm3']=all(v<=.001 for row in report['actual_frame_material_intersections_mm3'].values() for v in row.values())
    report['not_verified']=['three contact pads/jaws and filled wedge gaps','base mounting fasteners',
        'printability, stiffness and strength','assembly access and removable fastening',
        'continuous bone or deformed-muscle motion','GH retention and real fixture force path']
    (folder/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Lateral base trial count',len(trials),'selected',chosen[1]['displacement_source_mm'] if chosen else None,flush=True)
    if chosen:print('Mesh',{k:v['mesh_checks_pass'] for k,v in report['parts'].items()},
        'actual frame hits',report['actual_frame_material_intersections_mm3'],flush=True)


if __name__=='__main__':main()
