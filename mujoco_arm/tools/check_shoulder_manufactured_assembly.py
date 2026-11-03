"""Actual mesh assembly and local port probes. No blanket contact exclusions.

Contains available parts only; missing fixture jaws, screws and liner are
reported, not replaced by fictitious constraints or treated as validated.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import numpy as np
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0,str(Path(__file__).resolve().parent))
import frames
import interference as I
import mesh_geom
import search_shoulder_channel_paths as SP
import build_shoulder_tendon_segments as T

ROOT=frames.PROJECT_ROOT
DESIGN=frames.OUTPUT_DIR.parent/'shoulder_design'
REV=DESIGN/'s1_integrated_revision'


def load(path,translation=None):
    V,F=mesh_geom.read_stl(str(path));W,G=mesh_geom.weld(V,F)
    s=I._manifold().Manifold(I._manifold().Mesh64(W,G.astype(np.uint32)))
    if s.volume()<=0:raise ValueError('Invalid or empty mesh: '+str(path))
    return s if translation is None else s.translate(translation)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--unified-report',type=Path,action='append',help='Replace a complete module; repeat for multiple modules')
    ap.add_argument('--fixture-report',type=Path,help='Use checked actual fixture solids instead of the five archived pieces')
    ap.add_argument('--output-dir',type=Path,help='Write a separate candidate assembly report')
    args=ap.parse_args()
    report_paths=args.unified_report or []
    unified_reports=[(p,json.loads(p.read_text(encoding='utf8'))) for p in report_paths]
    replaced={r['module'] for _,r in unified_reports}
    if len(replaced)!=len(unified_reports):raise ValueError('Duplicate module replacement')
    dest=report_paths[0].parent if len(report_paths)==1 else DESIGN/'s1_port_path_revision/whole_shoulder_static' if report_paths else REV
    if args.fixture_report:dest=DESIGN/'s1_port_path_revision/whole_shoulder_fixture_revision'
    if args.output_dir:dest=args.output_dir
    dest.mkdir(exist_ok=True)
    parts={};sources={};kinds={};owners={}
    def add(key,file,kind,owner=None,translation=None):
        file=Path(file);s=SP.load_solid(file) if file.suffix.lower()=='.obj' else load(file,translation)
        parts[key]=s;kinds[key]=kind;owners[key]=owner
        sources[key]=dict(file=str(file.relative_to(ROOT)),sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                          translation_source_mm=translation,material_components=len([c for c in s.decompose() if c.volume()>1e-6]))
    bone_dir=ROOT/'02_人体骨骼晶格肌肉方案/输入/骨骼'
    for name,file in SP.BONES.items():add('bone:'+name,bone_dir/file,'bone')
    for name in ['deltoid','infraspinatus','subscapularis','teres_major']:
        if name in replaced:continue
        add('tpu:'+name,frames.OUTPUT_DIR/'shoulder_tpu'/(name+'_tpu.stl'),'TPU',name)
    seats=json.loads((DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    for s in seats:
        if s['module'] in replaced:continue
        filename=','.join(s['heads']).replace('/','_').replace(',','+')
        label=s['head']+'/'+s['end']
        for kind in ['pad','cap']:
            add(kind+':'+label,frames.OUTPUT_DIR/'shoulder_parts'/(filename+'_'+kind+'.stl'),kind,s['module'])
    segments=json.loads((REV/'tendons_export.json').read_text(encoding='utf8'))['segments']
    for row in segments:
        if row['built'] and row['module'] not in replaced:
            add('tendon:'+row['head']+'/'+row['end'],REV/'tendons'/row['file'],'tendon',row['module'],row['translation_source_mm'])
    fixture=json.loads(args.fixture_report.read_text(encoding='utf8')) if args.fixture_report else None
    if fixture:
        if not fixture['parts']:raise ValueError('Empty fixture replacement')
        for key,row in fixture['parts'].items():
            file=ROOT/row['file']
            if not row['mesh_checks_pass'] or hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError('Invalid or changed fixture replacement '+key)
            add('fixture:'+key,file,'fixture',translation=row['translation_source_mm'])
    else:
        for f in sorted((frames.OUTPUT_DIR/'shoulder_parts_fixture').glob('*.stl')):
            add('fixture:'+f.stem,f,'fixture')
    master_routes=[]
    master_hashes=set()
    for report_path,unified in unified_reports:
        module=unified['module']
        for key,row in unified['parts'].items():
            file=ROOT/row['file']
            if hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError('Unified part changed since its check: '+str(file))
            is_tpu=key=='TPU';is_stop=row.get('kind')=='stop_blank'
            label='tpu:'+module if is_tpu else ('stop:' if is_stop else 'cap:')+key
            add(label,file,'TPU' if is_tpu else 'stop_blank' if is_stop else 'cap',module,row['translation_source_mm'])
        source=ROOT/unified['master_routes']
        if hashlib.sha256(source.read_bytes()).hexdigest()!=unified['master_sha256']:
            raise ValueError('Unified master routes changed since generation')
        master_hashes.add(unified['master_sha256'])
        if len(master_hashes)>1:raise ValueError('Cannot combine different master route revisions')
        local_routes=[r for r in json.loads(source.read_text(encoding='utf8'))['routes'] if r['module']==module]
        master_routes.extend(local_routes)
        if module=='deltoid':
            folder=report_path.parent.parent/'triplex_terminal'
            terminal=json.loads((folder/'check.json').read_text(encoding='utf8'))
            ref_file=ROOT/terminal['routes_source']
            if hashlib.sha256(ref_file.read_bytes()).hexdigest()!=terminal['routes_sha256']:
                raise ValueError('Triplex reference routes changed')
            reference={r['head']:r for r in json.loads(ref_file.read_text(encoding='utf8'))['routes']}
            for r in local_routes:
                if r['ports']['distal']!=reference[r['head']]['ports']['distal']:
                    raise ValueError('Cannot reuse stop blank: distal port changed')
            for key,row in terminal['parts'].items():
                if not key.endswith('_stop_blank'):continue
                file=folder/row['file']
                if hashlib.sha256(file.read_bytes()).hexdigest()!=row['sha256']:
                    raise ValueError('Triplex stop blank changed')
                add('stop:'+key,file,'stop_blank',module,row['translation_source_mm'])
    contacts=[]
    nonzero_pairs=[]
    for a,b in itertools.combinations(parts,2):
        volume=I.intersection_volume(parts[a],parts[b])
        if volume>1e-12:nonzero_pairs.append(dict(a=a,b=b,volume_mm3=volume))
        if volume<=.05:continue
        same_owner=owners[a] is not None and owners[a]==owners[b]
        contacts.append(dict(a=a,b=b,volume_mm3=volume,
            classification='same-module connection requires fit verification' if same_owner else
                           'compliant-contact review' if 'TPU' in [kinds[a],kinds[b]] or 'pad' in [kinds[a],kinds[b]] else
                           'rigid overlap requires fit verification'))
    probes=[]
    for row in segments:
        if not row['built'] or row['module'] in replaced:continue
        P=np.array(row['path_source_mm']);rad=row['outer_radius_mm']
        a=(P[1]-P[0])/np.linalg.norm(P[1]-P[0]);b=(P[-1]-P[-2])/np.linalg.norm(P[-1]-P[-2])
        points=np.vstack([P[0]-a*(rad+1),P,P[-1]+b*(rad+1)])
        pieces=[T.rod(p,q,.3) for p,q in zip(points,points[1:])]
        pieces += [I._manifold().Manifold.sphere(.3,16).translate(p) for p in points]
        wire=I._manifold().Manifold.batch_boolean(pieces,I._manifold().OpType.Add)
        core_pieces=[T.rod(p,q,.3) for p,q in zip(P,P[1:])]
        core_pieces += [I._manifold().Manifold.sphere(.3,16).translate(p) for p in P]
        core=I._manifold().Manifold.batch_boolean(core_pieces,I._manifold().OpType.Add)
        own='tendon:'+row['head']+'/'+row['end']
        hits=[];core_hits=[]
        for name,s in parts.items():
            v=I.intersection_volume(wire,s)
            if v>.001:hits.append(dict(part=name,volume_mm3=v))
            vcore=I.intersection_volume(core,s)
            if vcore>.001:core_hits.append(dict(part=name,volume_mm3=vcore))
        probes.append(dict(head=row['head'],end=row['end'],own_tendon_overlap_mm3=I.intersection_volume(wire,parts[own]),
                           material_hits=hits,known_centreline_hits=core_hits,
                           scope='0.6mm nominal wire; separate known tendon centreline from exploratory straight port extensions'))
    master_probes=[]
    if master_routes:
        import design_shoulder_port_paths as D
        for row in master_routes:
            probe=D.outer_envelope(np.array(row['full_wire_path_source_mm']),.3)
            hits=[dict(part=k,volume_mm3=I.intersection_volume(probe,s)) for k,s in parts.items()]
            master_probes.append(dict(head=row['head'],all_material_intersections_mm3=hits,
                material_hits=[h for h in hits if h['volume_mm3']>.001],scope='complete nominal master steel path, all available parts'))
    report=dict(scope='actual available manufacturing meshes at zero pose, NOT full assembly or motion acceptance',
        n_instances=len(parts),sources=sources,contacts=contacts,wire_probes=probes,
        master_wire_probes=master_probes,replaced_modules=sorted(replaced),
        input_mesh_checks={r['module']:{k:v['mesh_checks_pass'] for k,v in r['parts'].items()} for _,r in unified_reports},
        fixture_input_sha256=hashlib.sha256(args.fixture_report.read_bytes()).hexdigest() if args.fixture_report else None,
        nonzero_material_pair_intersections_mm3=nonzero_pairs,
        unified_report_sha256={str(p.relative_to(ROOT) if p.is_absolute() else p):hashlib.sha256(p.read_bytes()).hexdigest() for p in report_paths},
        missing=['three fixture contact pads not exported','fixture jaws/bolts/rails','muscle fasteners/crimps','GH2mm liner','upstream conduit routing and reaction mounts'],
        not_verified=['attachment and retention','flexible deformation','wire minimum bend radius and friction','moving paths'],stage_pass=False)
    (dest/'actual_assembly_check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Actual instances',len(parts),'contacts >0.05mm3',len(contacts))
    print('Wire probes with other material hits',sum(bool(p['material_hits']) for p in probes),'/',len(probes))
    for p in probes:
        print(p['head'],p['end'],'own tube hit',round(p['own_tendon_overlap_mm3'],6),
              'other hits',[(r['part'],round(r['volume_mm3'],3)) for r in p['material_hits']])
        print('  known centreline hits',[(r['part'],round(r['volume_mm3'],3)) for r in p['known_centreline_hits']])
    for p in master_probes:print('Master',p['head'],'hits',p['material_hits'])


if __name__=='__main__':main()
