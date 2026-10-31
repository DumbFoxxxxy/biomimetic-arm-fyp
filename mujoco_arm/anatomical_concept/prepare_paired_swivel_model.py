"""Three actual hinge coordinates, rigid-only prescribed geometry demonstration.

No actuator, TPU deformation or contact physics is inferred from this viewer.
"""
from build_direct_muscle_yoke import *
import xml.etree.ElementTree as ET

def run():
    dest=DEST/'paired_swivel_support'
    if '--keepers' in sys.argv:dest=dest/'captive_keepers'
    a=json.loads((dest/'assembly.json').read_text());j=json.loads((BASE/'design.json').read_text());E=np.array(j['joint_centre_display_mm']);root=ET.Element('mujoco',model='paired_caps_rigid_prescribed_geometry_only');ET.SubElement(root,'compiler',angle='degree');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 2',dir='0 1 -2');forearm=ET.SubElement(world,'body',name='forearm',pos=' '.join(map(str,E/1000)));ET.SubElement(forearm,'joint',name='elbow',type='hinge',axis=' '.join(map(str,j['joint_axis'])),range='0 150',limited='true');owners={'elbow':(forearm,E)}
    for fid,folder in [('FJ1486',DEST/'cap_swivel_flexor'),('FJ1477',DEST/'cap_swivel_study')]:
        s=json.loads((folder/'screen.json').read_text());C=np.array(s['pivot_zero_display_mm']);body=ET.SubElement(forearm,'body',name=fid+'_distal_pivot',pos=' '.join(map(str,(C-E)/1000)));ET.SubElement(body,'joint',name=fid+'_swivel',type='hinge',axis=' '.join(map(str,s['trunnion_axis_zero_display'])),limited='false');owners[fid+'_distal_pivot']=(body,C)
    included=[];omitted=[]
    for r in a['instances']:
        if r['kind']=='TPU' or r.get('nominal_packing_only'):omitted.append(r['key']);continue
        key=r['key'];included.append(key);ET.SubElement(asset,'mesh',name=key,file=str(ROOT/r['file']),scale='.001 .001 .001');body,centre=owners.get(r['owner'],(world,np.zeros(3)));colour='.2 .55 .72 1' if key=='paired_bone_frame' else '.8 .57 .15 1' if key.endswith('_cap_rotor') else '.8 .81 .76 1';ET.SubElement(body,'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,(np.array(r['translation_display_mm'])-centre)/1000)),rgba=colour,contype='0',conaffinity='0',mass='.001')
    ET.indent(root);ET.ElementTree(root).write(dest/'rigid_hinges.xml',encoding='utf-8',xml_declaration=True);(dest/'rigid_model_scope.json').write_text(json.dumps({'scope':f'{len(included)}rigid parts, elbow plus2passive cap hinges; prescribed orientation demonstration only. Placeholder masses and disabled simulator contacts; actual solid geometry checked independently','included':included,'omitted':omitted,'joint_names':['elbow','FJ1486_swivel','FJ1477_swivel'],'drive_verified':False,'frozen':False},indent=2),encoding='utf8');print('paired hinge model',len(included),'rigid parts,omitted',len(omitted),flush=True)

if __name__=='__main__':run()
