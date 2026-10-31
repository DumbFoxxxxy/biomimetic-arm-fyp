"""Integrate candidate elbow geometry into full neutral A1, without overwrite."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json,xml.etree.ElementTree as ET
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_arm/A1'

def run():
    dest=OUT/'elbow_connection';tree=ET.parse(OUT/'arm.xml');root=tree.getroot();root.set('model','A1_connected_elbow_candidate');root.find('compiler').set('meshdir','.')
    for mesh in root.findall('./asset/mesh'):mesh.set('file','../assets/'+mesh.get('file'))
    for n,file in [('humerus','humerus_bored.stl'),('ulna','forearm_bridged.stl'),('elbow_pin','retained_axle.stl')]:root.find(f'./asset/mesh[@name="{n}"]').set('file',file)
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='geom' and child.get('name') in ['radius','elbow_bearing_reference']:parent.remove(child)
    frames=json.loads((OUT/'body_frames.json').read_text(encoding='utf8'))['origins_mm']
    for n,owner in [('cotter_pin','humerus'),('thrust_washer_-1','forearm'),('thrust_washer_1','forearm')]:
        ET.SubElement(root.find('asset'),'mesh',name=n,file=n+'.stl',scale='.001 .001 .001')
        ET.SubElement(root.find(f'.//body[@name="{owner}"]'),'geom',name=n,mesh=n,type='mesh',pos=' '.join(str(-v/1000) for v in frames[owner]),rgba='.25 .5 .65 1',contype='0',conaffinity='0',group='3',mass='.001')
    ET.indent(root);tree.write(dest/'arm.xml',encoding='utf-8',xml_declaration=True)
    print(dest/'arm.xml')

if __name__=='__main__':run()
