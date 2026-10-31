"""Independent 10-DOF joint-constraint prototype; no tendons/forces yet.
SC/AC/GH spherical joints are mathematical constraints, not invented CAD caps.
"""
from build_reference import *

def run():
    fit=json.loads((OUT/'joint_candidates.json').read_text(encoding='utf8'));inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));dest=OUT/'joint_model';dest.mkdir(exist_ok=True)
    centres={k:np.array(fit[k]['centre_source_mm']) for k in ['SC','AC','GH','elbow']};root=ET.Element('mujoco',model='independent_SC3_AC3_GH3_elbow1_candidate');ET.SubElement(root,'compiler',angle='degree',meshdir='../assets');ET.SubElement(root,'option',gravity='0 0 0');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='1000');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1');nodes={};owners={'FJ3362':'SC','FJ3384':'AC','FJ3368':'GH','FJ3391':'elbow','FJ3349':'elbow'}
    for name,parent in [('SC',None),('AC','SC'),('GH','AC'),('elbow','GH')]:
        offset=centres[name]-(centres[parent] if parent else SHIFT);node=ET.SubElement(nodes[parent] if parent else world,'body',name=name,pos=' '.join(str(x/1000) for x in offset));nodes[name]=node
        if name=='elbow':ET.SubElement(node,'joint',name='elbow_hinge',type='hinge',axis=' '.join(map(str,fit[name]['axis_source'])),range='-5 170',damping='1')
        else:ET.SubElement(node,'joint',name=name+'_spherical',type='ball',damping='1')
        ET.SubElement(node,'site',name=name+'_centre',pos='0 0 0',size='.003',rgba='.2 .5 1 1',group='3')
    for rec in inv['models']:
        fid=rec['id'];ET.SubElement(asset,'mesh',name=fid,file=fid+'.stl',scale='.001 .001 .001');attrs={'name':fid,'mesh':fid,'type':'mesh','contype':'0','conaffinity':'0'}
        if fid in owners:
            owner=owners[fid];attrs.update(pos=' '.join(str(x/1000) for x in SHIFT-centres[owner]),rgba='.9 .87 .8 1',mass='.1',group='0');ET.SubElement(nodes[owner],'geom',**attrs)
        else:
            attrs.update(rgba='.75 .12 .12 .35' if rec['kind']=='muscle' else '.7 .75 .8 .3',group='1' if rec['kind']=='muscle' else '2');ET.SubElement(world,'geom',**attrs)
    ET.indent(root);ET.ElementTree(root).write(dest/'arm.xml',encoding='utf8',xml_declaration=True)
    report={'scope':'engineering constraint candidate from independently fitted raw source anatomy; not manufactured or driven assembly','joints':{'SC':'spherical3','AC':'spherical3','GH':'spherical3','elbow':'hinge1'},'velocity_dofs':10,'frozen':False,'centres_source_mm':{k:v.tolist() for k,v in centres.items()},'bone_owner':owners,'source_coordinate_shift_mm':SHIFT.tolist(),'muscles':'neutral reference surfaces only; must hide while changing joints until compliant routing model is defined','forearm':'source ulna/radius held together solely for first elbow geometric test; final pronation/supination not determined','not_completed':['physical joint capture and supports','ligament/contact restraint of scapula','muscle motion and steel channels','force/strength/printed performance']};(dest/'definition.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print(dest/'arm.xml',flush=True)

if __name__=='__main__':run()
