"""Independent raw-coordinate anatomy reference; no inherited joint chain.
All anatomy meshes undergo the same translation and no scaling or mirroring.
"""
from pathlib import Path
import sys,json,hashlib,struct,re,itertools
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'骨骼/_optimize'))
from decimate import read_obj,weld_stl
import xml.etree.ElementTree as ET
OUT=ROOT/'mujoco_arm/exports/anatomical_concept'
SHIFT=np.array([-140.,-70.,1100.])

def run():
    (OUT/'assets').mkdir(parents=True,exist_ok=True);rows=[]
    root=ET.Element('mujoco',model='independent_anatomical_reference');ET.SubElement(root,'compiler',meshdir='assets');ET.SubElement(root,'option',gravity='0 0 0');visual=ET.SubElement(root,'visual');ET.SubElement(visual,'global',offwidth='1200',offheight='1000');asset=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 1',dir='0 1 -1')
    paths=list((ROOT/'骨骼').glob('*.obj'))+list((ROOT/'muscle').glob('*.obj'))+list((ROOT/'骨骼/支撑参考').glob('*.obj'))
    for path in paths:
        header=path.read_text(encoding='utf8',errors='replace').splitlines()[:25];metadata={}
        for line in header:
            if line.startswith('# ') and ' : ' in line:
                k,v=line[2:].split(' : ',1);metadata[k.strip()]=v.strip()
        fid=metadata.get('File ID',path.stem);name=metadata.get('English name',path.stem);v,f=read_obj(str(path));rawbounds=[v.min(0).tolist(),v.max(0).tolist()];v,f=weld_stl(v[f]);edges=Counter(tuple(sorted((int(a),int(b)))) for face in f for a,b in [(face[0],face[1]),(face[1],face[2]),(face[2],face[0])]);boundary=sum(n==1 for n in edges.values());nonmanifold=sum(n>2 for n in edges.values());v=v-SHIFT
        tri=v[f];norm=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm/=np.maximum(np.linalg.norm(norm,axis=1),1e-30)[:,None];dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);data=np.zeros(len(f),dt);data['n']=norm;data['v']=tri
        output=OUT/'assets'/(fid+'.stl');output.write_bytes(b'Unmodified source anatomy, shared translation only, mm'.ljust(80,b' ')+struct.pack('<I',len(f))+data.tobytes());kind='bone' if path.parent.name=='骨骼' else 'muscle' if path.parent.name=='muscle' else 'support_reference'
        ET.SubElement(asset,'mesh',name=fid,file=fid+'.stl',scale='.001 .001 .001')
        ET.SubElement(world,'geom',name=fid,mesh=fid,type='mesh',rgba='.9 .87 .8 1' if kind=='bone' else '.75 .12 .12 .6' if kind=='muscle' else '.7 .75 .8 .35',group='0' if kind=='bone' else '1' if kind=='muscle' else '2',contype='0',conaffinity='0')
        rows.append({'id':fid,'name':name,'kind':kind,'source':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'metadata':metadata,'raw_bounds_mm':rawbounds,'vertices':len(v),'triangles':len(f),'boundary_edges':boundary,'nonmanifold_edges':nonmanifold,'mesh_closed_edge_test':boundary==0 and nonmanifold==0,'attachment_regions':'not yet labeled','fiber_directions':'not supplied','geometry_status':'anatomical surface only, not TPU lattice'})
    ET.indent(root);ET.ElementTree(root).write(OUT/'reference.xml',encoding='utf8',xml_declaration=True)
    required=['deltoid','pectoralis major','trapezius','serratus anterior','rhomboid'];names=' '.join(r['name'].lower() for r in rows if r['kind']=='muscle')
    report={'route':'independent anatomical concept; no V2/A1 kinematic constraints','translation_subtracted_mm':SHIFT.tolist(),'rotation':'identity','scale':1,'laterality':'Right, from source labels','units':'mm, from source header','registration':'common raw BodyParts3D coordinates retained; header compatibility alone does not prove physiological fit','models':rows,'license':'BodyParts3D, (c) The Database Center for Life Science, CC Attribution-Share Alike 2.1 Japan, as stated in source OBJ headers','source_license_url':'http://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html','unprovided_or_unidentified_muscles':[name for name in required if name not in names],'required_group_presence_not_full_anatomical_completeness':True,'drive_validation':False}
    (OUT/'inventory.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    txt=['# 人体骨肌输入清单','','本模型保留所有原始公共坐标，整体平移至方便查看的位置，未拟合或缩放单块肌肉。','这是解剖外形参考，不是晶格或驱动仿真。','','|ID|源文件英文名称|类型|三角面|闭合边检查|','|---|---|---|---:|---|']
    for row in rows:txt.append(f"|{row['id']}|{row['name']}|{row['kind']}|{row['triangles']}|{row['mesh_closed_edge_test']}|")
    txt+=['','Source attribution: BodyParts3D, (c) The Database Center for Life Science,','CC Attribution-Share Alike 2.1 Japan. Source header license URL is preserved in inventory.json.','','现有模型没有标注肌纤维方向及精确附着区域。缺少或尚未识别的关键肌群见inventory.json。','来源名称已确认是右侧；坐标叠合仍需视觉及局部接触核验。']
    (OUT/'INVENTORY.md').write_text('\n'.join(txt)+'\n',encoding='utf8');print('bones',sum(r['kind']=='bone' for r in rows),'muscles',sum(r['kind']=='muscle' for r in rows),'nonclosed',[(r['id'],r['boundary_edges'],r['nonmanifold_edges']) for r in rows if not r['mesh_closed_edge_test']],flush=True)

if __name__=='__main__':run()
