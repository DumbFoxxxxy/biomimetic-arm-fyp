"""Independent M7 neutral and determinate rigid models; no cap-motion law."""
from build_guide_ground_frames import *
import os

def run():
    design=json.loads((DEST7/'design.json').read_text());assert all(m['selected'] for m in design['modules'])
    a=json.loads((CURRENT/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};scope=json.loads((BASE/'rigid_motion_scope.json').read_text());tree=ET.parse(BASE/'arm.xml');root=tree.getroot();root.set('model','M7_known_rigid_guide_frames_only')
    root.find('compiler').set('meshdir','.')
    for mesh in root.findall('./asset/mesh'):mesh.set('file',str((BASE/mesh.get('file')).resolve()))
    free_stops=['flexor_tendon_stop','extensor_tendon_stop']
    for name in free_stops:
        rows[name]['owner']='distal_unvalidated'
        for parent in root.iter():
            for child in list(parent):
                if child.tag in ['geom','mesh'] and child.get('name')==name:parent.remove(child)
    gh=np.array(json.loads((OUT/'joint_model/definition.json').read_text())['centres_source_mm']['GH'])-np.array(json.loads((OUT/'joint_model/definition.json').read_text())['source_coordinate_shift_mm'])
    manufactured_names={p['name'] for p in json.loads((DEST7/'cad_input.json').read_text())['parts']}
    for fid in ['FJ1486','FJ1477']:
        key=fid+'_bone_adapter';rows[key]['file']=str((DEST7/(key+'.obj')).relative_to(ROOT));rows[key]['translation_display_mm']=[0,0,0];root.find(f'./asset/mesh[@name="{key}"]').set('file',str(DEST7/(key+'.obj')))
        for path in sorted((DEST/'paired_placement').glob(fid+'_*.obj')):
            if path.stem.endswith('wire_corridor'):continue
            rows[path.stem]={'key':path.stem,'file':str(path.relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':'GH','kind':'rigid','grounded_guide_candidate':True}
        for path in [DEST7/(name+'.obj') for name in sorted(manufactured_names) if name.startswith(fid+'_key_')]:
            rows[path.stem]={'key':path.stem,'file':str(path.relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':'GH','kind':'rigid','grounded_guide_candidate':True}
    existing={g.get('name') for g in root.findall('.//geom')};included=[];omitted=[]
    for key,r in rows.items():
        if r['owner'] in ['soft_unvalidated','distal_unvalidated']:omitted.append(key);continue
        included.append(key)
        if key in existing:continue
        assert r['owner']=='GH'
        ET.SubElement(root.find('asset'),'mesh',name=key,file=str(ROOT/r['file']),scale='.001 .001 .001')
        ET.SubElement(root.find('.//body[@name="GH"]'),'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,-gh/1000)),rgba='.2 .5 .8 1',group='3',contype='0',conaffinity='0',mass='.001')
    a.update(scope='M7 neutral manufactured pieces including unverified free muscle ends; only grounded frame/guide pieces enter prescribed rigid motion model',instances=list(rows.values()),frozen=False)
    (DEST7/'assembly.json').write_text(json.dumps(a,indent=2),encoding='utf8');ET.indent(root);tree.write(DEST7/'arm.xml',encoding='utf-8',xml_declaration=True)
    (DEST7/'rigid_motion_scope.json').write_text(json.dumps({'scope':'M7 known rigid bones, nominal bolts, integral adapter extensions and key-mounted sheaves; no installed cap/TPU/transmission pass','included':included,'omitted_compliant_or_free_distal':omitted,'omitted_old_uncrimped_stops':free_stops,'frozen':False},indent=2),encoding='utf8')
    print('M7 neutral',len(rows),'determinate rigid',len(included),'omitted',len(omitted),flush=True)

if __name__=='__main__':run()
