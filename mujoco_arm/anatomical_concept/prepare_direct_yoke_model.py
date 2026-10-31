"""Direct cap-to-bone boundaries are actual bolted geometry; TPU unsolved."""
from build_direct_muscle_yoke import *
import xml.etree.ElementTree as ET

def run():
    a=json.loads((DEST/'assembly.json').read_text());j=json.loads((OUT/'joint_model/definition.json').read_text());centres={k:np.array(v)-np.array(j['source_coordinate_shift_mm']) for k,v in j['centres_source_mm'].items()};tree=ET.parse(M2/'arm.xml');root=tree.getroot();root.set('model','M9_direct_TPU_bone_interfaces_rigid_only');root.find('compiler').set('meshdir','.')
    for asset in root.findall('./asset/mesh'):asset.set('file',str((M2/asset.get('file')).resolve()))
    included=[];omitted=[];boundaries=[]
    for r in a['instances']:
        if r.get('nominal_packing_only'):r['nominal_pose_owner']='elbow';r['owner']='distal_unvalidated'
        key=r['key']
        if r['owner'] in ['soft_unvalidated','distal_unvalidated']:omitted.append(key);continue
        included.append(key);asset=root.find(f'./asset/mesh[@name="{key}"]')
        if asset is None:asset=ET.SubElement(root.find('asset'),'mesh',name=key,scale='.001 .001 .001')
        asset.set('file',str(ROOT/r['file']));geom=root.find(f'.//geom[@name="{key}"]');parent=root.find('worldbody') if r['owner']=='world' else root.find(f'.//body[@name="{r["owner"]}"]')
        if geom is None:geom=ET.SubElement(parent,'geom',name=key,type='mesh',mesh=key,rgba='.2 .55 .7 1',group='3',contype='0',conaffinity='0',mass='.001')
        c=np.zeros(3) if r['owner']=='world' else centres[r['owner']];geom.set('pos',' '.join(map(str,(np.array(r['translation_display_mm'])-c)/1000)))
    for fid in ['FJ1486','FJ1477']:
        folder=muscle_folder(fid);d=json.loads((folder/'design.json').read_text());offset=np.array(next(p for p in a['placements'] if p['id']==fid)['translation_display_mm']);boundaries.append({'muscle':fid,'proximal_owner':'GH','distal_owner':'elbow','proximal_pad_centre_display_mm':(np.array(d['pad_centres_display_mm'][0])+offset).tolist(),'distal_pad_zero_centre_display_mm':(np.array(d['pad_centres_display_mm'][1])+offset).tolist(),'orientation_basis_columns_at_zero':d['basis_columns_display'],'constraint_basis':'actual shared M3 caps and one-piece bone yoke; clamp/preload/TPU slip not verified','soft_mesh_deformation_solved':False})
    a['muscle_end_boundaries']=boundaries;a['pending']=['rigid motion and source-to-CAD checks','TPU nonlinear deformation, bone contact,guide bending andcrimp grip','extensor joint-covering path andactual force sign','printing/material calibration','whole shoulder/scapular actuation andmotor reaction frame'];(DEST/'assembly.json').write_text(json.dumps(a,indent=2),encoding='utf8');ET.indent(root);tree.write(DEST/'arm.xml',encoding='utf8',xml_declaration=True);(DEST/'rigid_motion_scope.json').write_text(json.dumps({'scope':'direct bone-cap rigid bodies only, actual muscle-end mechanical boundaries; no TPU deformation law','included':included,'omitted':omitted,'uncrimped_nominal_stops':a['uncrimped_nominal_stops'],'frozen':False},indent=2),encoding='utf8');print('M9 direct rigid',len(included),'neutral',len(a['instances']),'omitted',len(omitted),flush=True)

if __name__=='__main__':run()
