"""Integrate one M4 candidate as111 neutral/89 known-rigid instances."""
from build_elbow_side_lever import *
import os,xml.etree.ElementTree as ET
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15';folder=DEST/variant
    design=json.loads((folder/'design.json').read_text(encoding='utf8'));assembly=json.loads((MOUNT/'assembly.json').read_text(encoding='utf8'))
    for r in assembly['instances']:
        if r['key']=='FJ3391':r['file']=str((folder/'forearm_spotfaced.obj').relative_to(ROOT));r['translation_display_mm']=[0,0,0]
    for name in design['parts']:
        if name=='forearm_spotfaced':continue
        assembly['instances'].append({'key':name,'file':str((folder/(name+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':'elbow','kind':'rigid','generated':True})
    assembly['generated_components']=design['parts'];assembly['scope']='M4 removable side lever candidate; no solved tendon or TPU deformation'
    (folder/'assembly.json').write_text(json.dumps(assembly,indent=2),encoding='utf8')
    tree=ET.parse(MOUNT/'arm.xml');root=tree.getroot();root.set('model','M4_known_rigid_side_lever_only')
    for mesh in root.findall('./asset/mesh'):mesh.set('file',os.path.relpath((MOUNT/mesh.get('file')).resolve(),folder).replace('\\','/'))
    root.find('./asset/mesh[@name="FJ3391"]').set('file','forearm_spotfaced.obj')
    E=np.array(design['joint_centre_display_mm'])
    for name in design['parts']:
        if name=='forearm_spotfaced':continue
        ET.SubElement(root.find('asset'),'mesh',name=name,file=name+'.obj',scale='.001 .001 .001')
        ET.SubElement(root.find('.//body[@name="elbow"]'),'geom',name=name,type='mesh',mesh=name,pos=' '.join(str(c/1000) for c in -E),rgba='.15 .55 .8 1' if name=='side_lever' else '.5 .6 .65 1',group='3',contype='0',conaffinity='0',mass='.001')
    ET.indent(root);tree.write(folder/'arm.xml',encoding='utf8',xml_declaration=True)
    scope=json.loads((MOUNT/'rigid_motion_scope.json').read_text(encoding='utf8'));scope['included']+=[k for k in design['parts'] if k!='forearm_spotfaced'];scope['scope']='M4 prescribed rigid packing:87 nominal bone/bolted pieces plus2 uncrimped passive stop positions; twoTPU and20 free distal pieces omitted, no grounding or drive proof for stops'
    scope['prescribed_uncrimped_stops']=['flexor_tendon_stop','extensor_tendon_stop'];scope['nominally_bolted_or_bone_instances']=87
    (folder/'rigid_motion_scope.json').write_text(json.dumps(scope,indent=2),encoding='utf8');print('M4 neutral',len(assembly['instances']),'knownrigid',len(scope['included']),flush=True)
if __name__=='__main__':run()
