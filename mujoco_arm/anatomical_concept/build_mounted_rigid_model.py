"""M2 known rigid grounding motion model; unvalidated compliant ends omitted."""
from build_elbow_muscle_mounts import *

def run():
    data=json.loads((DEST/'assembly.json').read_text());definition=json.loads((OUT/'joint_model/definition.json').read_text());G=np.array(definition['centres_source_mm']['GH'])-SHIFT
    tree=ET.parse(SKELETON/'arm.xml');root=tree.getroot();root.set('model','M2_known_rigid_proximal_grounding_only');root.find('compiler').set('meshdir','.')
    for mesh in root.findall('./asset/mesh'):mesh.set('file',os.path.relpath((SKELETON/mesh.get('file')).resolve(),DEST).replace('\\','/'))
    hum=next(r for r in data['instances'] if r['key']=='FJ3368');root.find('./asset/mesh[@name="FJ3368"]').set('file',os.path.relpath(ROOT/hum['file'],DEST).replace('\\','/'))
    existing={g.get('name') for g in root.findall('.//geom')};included=[];omitted=[]
    for r in data['instances']:
        if r['owner'] in ['soft_unvalidated','distal_unvalidated']:
            omitted.append(r['key']);continue
        included.append(r['key'])
        if r['key'] in existing:continue
        if r['owner']!='GH':raise RuntimeError('Unexpected new owner '+r['key'])
        key=r['key'];ET.SubElement(root.find('asset'),'mesh',name=key,file=os.path.relpath(ROOT/r['file'],DEST).replace('\\','/'),scale='.001 .001 .001')
        position=(np.array(r['translation_display_mm'])-G)/1000
        ET.SubElement(root.find('.//body[@name="GH"]'),'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,position)),rgba='.2 .45 .7 1' if 'adapter' in key or 'rigid_cap' in key else '.5 .6 .65 1',group='3',contype='0',conaffinity='0',mass='.001')
    ET.indent(root);tree.write(DEST/'arm.xml',encoding='utf-8',xml_declaration=True)
    (DEST/'rigid_motion_scope.json').write_text(json.dumps({'scope':'known rigid proximal parts only; not installed whole-module motion or muscle-drive acceptance','included':included,'omitted_compliant_or_free_distal':omitted,'omission_reason':'no validated TPU deformation/passive output transmission yet; no result transferred to these parts','frozen':False},indent=2),encoding='utf8');print('M2 rigid motion included',len(included),'omitted',len(omitted))

if __name__=='__main__':run()
