"""Penetration diagnostic of drawn tendon/cap paths; no rope force constraint.

Axial cap motion and bends at outlet references are UNVALIDATED assumptions.
This finds obstructions to proposed routes, never establishes flexible drive.
"""
from build_elbow_side_lever import *
from build_shaped_elbow_muscles import capsule
from scipy.spatial.transform import Rotation
def run():
    folder=DEST/'axis_negative_r15';design=json.loads((folder/'design.json').read_text(encoding='utf8'));assembly=json.loads((folder/'assembly.json').read_text(encoding='utf8'));diag=json.loads((folder/'tendon_axial_diagnostic.json').read_text(encoding='utf8'))
    rows={r['key']:r for r in assembly['instances']};local={r['key']:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for r in assembly['instances']}
    modules={fid:json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8')) for fid in ['FJ1486','FJ1477']};offsets={p['id']:np.array(p['translation_display_mm']) for p in assembly['placements']}
    paths={r['muscle']:r for r in diag['modules'] if (r['muscle'],r['lever_eyelet']) in [('FJ1486','flexor'),('FJ1477','extensor')]};assert len(paths)==2
    assert all(len(r['path'])==151 for r in paths.values())
    E=np.array(design['joint_centre_display_mm']);axis=np.array(design['joint_axis']);records=[]
    def overlap(a,b):
        aa=np.array(a.bounding_box());bb=np.array(b.bounding_box())
        if np.any(np.maximum(aa[:3]-bb[3:],bb[:3]-aa[3:])>=0):return 0.
        return max(0,float((a^b).volume()))
    for angle in range(151):
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();T=np.c_[R,E-R@E];world={};wires={};cap_groups={}
        for key,s in local.items():
            r=rows[key]
            if r['owner']=='elbow':s=s.transform(T)
            if r['owner']=='distal_unvalidated':
                fid=r['source_module'];z=np.array(modules[fid]['basis_columns_display'])[:,2];s=s.translate(-z*paths[fid]['path'][angle]['cap_axial_shortening_mm']);cap_groups.setdefault(fid,[]).append(key)
            world[key]=s
        row={'angle_deg':angle,'tendons':[],'cap_rigid_interferences':[]}
        for fid,path in paths.items():
            point=path['path'][angle];module=modules[fid];F=np.array(module['basis_columns_display']);short=point['cap_axial_shortening_mm'];pad=np.array(module['pad_centres_display_mm'][1])+offsets[fid]
            centre=pad+F@np.array([0,4,4.15])-F[:,2]*short;outlet=np.array(point['cap_exit_display_mm']);entry=np.array(point['lever_anchor_display_mm']);anchor0=np.array(design['eyelets_display_mm'][path['lever_eyelet']]);anchor=E+R@(anchor0-E);far=anchor+axis*design['axis_side']*5.55
            nodes=[centre,outlet,entry,far];wire=mf.Manifold.batch_boolean([capsule(a,b,.3) for a,b in zip(nodes,nodes[1:])],mf.OpType.Add);wires[fid]=wire;hard=[];soft=[]
            for key,s in world.items():
                volume=overlap(wire,s)
                if volume>.001:(soft if rows[key]['kind']=='TPU' else hard).append({'part':key,'volume_mm3':volume})
            row['tendons'].append({'muscle':fid,'prescribed_cap_shortening_mm':short,'nodes_display_mm':[p.tolist() for p in nodes],'rigid_path_penetrations':hard,'static_TPU_contact_flags':soft})
            for key in cap_groups[fid]:
                for other,s in world.items():
                    if other in cap_groups[fid] or rows[other]['kind']=='TPU':continue
                    volume=overlap(world[key],s)
                    if volume>.001:row['cap_rigid_interferences'].append({'pair':[key,other],'volume_mm3':volume})
        row['tendon_tendon_intersection_mm3']=overlap(wires['FJ1486'],wires['FJ1477']);records.append(row)
        if angle%30==0:print('tendon family',angle,'hard',sum(len(r['rigid_path_penetrations']) for r in row['tendons']),'cap',len(row['cap_rigid_interferences']),flush=True)
    summary=[]
    for fid in paths:
        hits=[{'angle_deg':r['angle_deg'],'hits':t['rigid_path_penetrations']} for r in records for t in r['tendons'] if t['muscle']==fid and t['rigid_path_penetrations']]
        summary.append({'muscle':fid,'rigid_penetration_samples':len(hits),'first_rigid_penetration':hits[0] if hits else None})
    report={'scope':'drawn tendon polyline and assumed axial rigid cap travel only; outlet-reference bends have no validated physical guiding/termination solution. Not an imposed MuJoCo constraint, deformation or drive acceptance.','step_deg':1,'wire_radius_mm':.3,'rigid_penetration_threshold_mm3':.001,'muscles':summary,'cap_interference_samples':sum(bool(r['cap_rigid_interferences']) for r in records),'wire_to_wire_interference_samples':sum(r['tendon_tendon_intersection_mm3']>.001 for r in records),'samples':records,'source_design_sha256':hashlib.sha256((folder/'design.json').read_bytes()).hexdigest(),'drive_pass':False,'frozen':False}
    (folder/'tendon_family_contacts.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(summary,'cap samples',report['cap_interference_samples'],flush=True)
if __name__=='__main__':run()
