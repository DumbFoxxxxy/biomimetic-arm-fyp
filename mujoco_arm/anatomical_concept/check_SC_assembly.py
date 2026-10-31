"""Validate actual exported SC geometry and quantify effects on candidate ROIs."""
from build_SC_connection import *
from extract_attachment_regions import projector
WITH_AC='--AC' in sys.argv
if WITH_AC:DEST=DEST/('AC_compact' if '--compact' in sys.argv else 'AC_connection')
WITH_GH='--GH' in sys.argv
if WITH_GH:DEST=OUT/'joint_model/SC_connection/AC_compact'/('GH_closed' if '--closed' in sys.argv else 'GH_preserve_head' if '--preserve-head' in sys.argv else 'GH_connection')
WITH_ELBOW='--elbow' in sys.argv
if WITH_ELBOW:DEST=OUT/'joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'

def run():
    data=json.loads((DEST/'cad_input.json').read_text(encoding='utf8'));checks=[];solids={}
    for part in data['parts']:
        name=part['name'];v,f=read_obj(str(DEST/(name+'.obj')));s=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)));solids[name]=s;raw=(DEST/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);re=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));err=abs(re.volume()-s.volume())/max(1,abs(s.volume()));coord=float(np.max(abs(tri-np.array(part['vertices'])[part['faces']])));checks.append({'part':name,'OBJ_status':str(s.status()),'STL_status':str(re.status()),'relative_volume_error':err,'STL_coordinate_error_mm':coord,'pass':s.status()==mf.Error.NoError and re.status()==mf.Error.NoError and err<1e-4 and coord<1e-4})
    collisions=[]
    for (an,aa),(bn,bb) in itertools.combinations(solids.items(),2):
        vol=max(0,float((aa^bb).volume()))
        if vol>.001:collisions.append({'parts':[an,bn],'volume_mm3':vol})
    targets={'clavicle':'clavicle_SC_AC_balls','scapula':'scapula_AC_socket'} if WITH_AC else {'clavicle':'clavicle_SC_ball'}
    if WITH_GH:targets={'scapula':'scapula_AC_GH_sockets','humerus':'humerus_GH_ball'}
    if WITH_ELBOW:targets={'humerus':'humerus_GH_elbow_fork','ulna':'forearm_elbow_hub','radius':'forearm_elbow_hub'}
    queries={}
    for structure,name in targets.items():
        v,f=read_obj(str(DEST/(name+'.obj')));queries[structure]=projector(v+SHIFT,f)
    atlasfile=OUT/'attachment_regions/atlas_and_regions.json';atlasdata=json.loads(atlasfile.read_text(encoding='utf8'));impact=[]
    for muscle in atlasdata['muscles']:
        for region in muscle['regions']:
            if region['structure'] not in queries or not region.get('selected_face_count_2mm'):continue
            pts=np.load(OUT/'attachment_regions'/region['data'])['bone_points_source_mm'];_,d,_=queries[region['structure']](pts);impact.append({'muscle':muscle['id'],'structure':region['structure'],'role':region['role'],'review_samples':len(d),'minimum_deviation_mm':float(d.min()),'median_deviation_mm':float(np.median(d)),'maximum_deviation_mm':float(d.max()),'points_more_than_1mm_from_new_surface':int(sum(d>1)),'frozen_anchor':False})
    report={'geometry_checks':checks,'neutral_all_pair_overlaps':collisions,'all_geometry_checks_pass':all(c['pass'] for c in checks),'old_anatomical_patch_impact':impact,'interpretation':'source ROI samples are not frozen anchors; any moved region needs remapping to actual engineered bone. Not a load/drive validation.'}
    if not WITH_AC and not WITH_GH and not WITH_ELBOW:report['old_anatomical_clavicle_patch_impact']=impact
    (DEST/'geometry_checks.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8');print('export geometry',len(checks),all(c['pass'] for c in checks),'neutral overlaps',collisions);print('patch impact',[(r['muscle'],r['structure'],r['points_more_than_1mm_from_new_surface'],r['maximum_deviation_mm']) for r in impact])

if __name__=='__main__':run()
