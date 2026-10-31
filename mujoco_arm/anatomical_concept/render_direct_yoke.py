"""Actual direct-muscle neutral assembly; steel is threading reference only."""
from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
import numpy as np,mujoco
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'mujoco_arm/exports/anatomical_concept';SHIFTED='--shifted' in sys.argv;DEST=OUT/'direct_muscle_interfaces'/('mechanical_yoke_shifted' if SHIFTED else 'mechanical_yoke_zero')
if '--access-relief' in sys.argv:DEST=OUT/'direct_muscle_interfaces/yoke_access_relief_study'
if '--extensor-reanchor' in sys.argv:DEST=OUT/'direct_muscle_interfaces/mechanical_yoke_reanchored_pair'

def run():
    a=json.loads((DEST/'assembly.json').read_text());root=ET.Element('mujoco',model='M9 direct muscle end interfaces, neutral only');ET.SubElement(root,'compiler',angle='radian');vis=ET.SubElement(root,'visual');ET.SubElement(vis,'global',offwidth='1400',offheight='900');assets=ET.SubElement(root,'asset');world=ET.SubElement(root,'worldbody');ET.SubElement(world,'light',pos='0 -1 2',dir='0 1 -2',diffuse='.85 .85 .85');hashes={}
    for r in a['instances']:
        key=r['key'];p=(ROOT/r['file']).resolve();hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();ET.SubElement(assets,'mesh',name=key,file=str(p),scale='.001 .001 .001');colour='.8 .81 .76 1' if key.startswith('FJ') else '.6 .65 .7 1'
        if r['kind']=='TPU':colour='.75 .1 .06 1'
        if key=='direct_muscle_yoke':colour='.2 .55 .72 1'
        ET.SubElement(world,'geom',name=key,type='mesh',mesh=key,pos=' '.join(map(str,np.array(r['translation_display_mm'])/1000)),rgba=colour,contype='0',conaffinity='0')
    segments=0
    for fid in ['FJ1486','FJ1477']:
        folder=OUT/'direct_muscle_interfaces/shifted_muscle'/fid if SHIFTED and fid=='FJ1486' else OUT/'elbow_muscles'/fid
        if '--extensor-reanchor' in sys.argv and fid=='FJ1477':folder=OUT/'direct_muscle_interfaces/extensor_reanchor_extended/muscle/FJ1477'
        d=json.loads((folder/'design.json').read_text());offset=np.array(next(p for p in a['placements'] if p['id']==fid)['translation_display_mm']);nodes=np.array(d['wire_polyline_display_mm'])+offset
        # Trim assembly feed overrun to the actual uncrimped stop midplane.
        z=np.array(d['basis_columns_display'])[:,2];nodes[-1]=np.array(d['pad_centres_display_mm'][1])+offset+z*4.6
        for i,(p,q) in enumerate(zip(nodes,nodes[1:])):ET.SubElement(world,'geom',name=f'wire_reference_{fid}_{i}',type='capsule',fromto=' '.join(map(str,np.r_[p,q]/1000)),size='.0003',rgba='1 .5 .05 1',contype='0',conaffinity='0');segments+=1
    xml=DEST/'neutral_view.xml';ET.ElementTree(root).write(xml,encoding='utf8',xml_declaration=True);m=mujoco.MjModel.from_xml_path(str(xml));d=mujoco.MjData(m);mujoco.mj_forward(m,d);renderer=mujoco.Renderer(m,height=900,width=1400)
    for name,lookat,distance in [('whole_arm',[-.05,0,.025],1.05),('direct_interface_detail',[-.07,0,-.07],.35)]:
        cam=mujoco.MjvCamera();cam.lookat[:]=lookat;cam.distance=distance;cam.azimuth=125;cam.elevation=-7;renderer.update_scene(d,camera=cam);renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW]=False;im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle((0,0,1400,64),fill='white');draw.text((12,10),'M9 | direct TPU end-to-forearm interfaces | 107 actual manufactured meshes | MuJoCo',fill='black');draw.text((12,34),'NEUTRAL only; internal steel reference; crimp/TPU deformation/actuation UNVERIFIED',fill='black');im.save(DEST/(name+'.png'))
    renderer.close();assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items());(DEST/'render_inputs.json').write_text(json.dumps({'scope':'107 actual neutral meshes plusinternal steel reference capsules, no drive/TPU deformation','meshes':len(a['instances']),'analytic_steel_reference_segments':segments,'source_sha256':hashes,'assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest(),'XML_sha256':hashlib.sha256(xml.read_bytes()).hexdigest()},indent=2),encoding='utf8');print('direct M9 render',len(a['instances']),'meshes +',segments,'wire reference segments',flush=True)

if __name__=='__main__':run()
