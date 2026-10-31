"""Separate rigid-housing print placements; metal blanks are not TPU prints."""
from build_distal_swivels import *
def run():
    results=[]
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=json.loads((folder/'design.json').read_text(encoding='utf8'));module=json.loads((OUT/'elbow_muscles'/fid/'design.json').read_text(encoding='utf8'));F=np.array(module['basis_columns_display']);output=folder/'print_candidates';output.mkdir(exist_ok=True)
        for name in ['distal_socket_base','socket_retainer']:
            z=F[:,2] if name=='distal_socket_base' else np.array(design['selected']['neutral_terminal_axis']);x=F[:,0] if name=='distal_socket_base' else np.array(design['selected']['retainer_v_axis']);x=unit(x-z*(x@z));y=np.cross(z,x);frame=np.stack([x,y,z],1)
            source,v,f=load(folder/(name+'.obj'));local=v@frame;shift=np.array([-(local[:,0].min()+local[:,0].max())/2,-(local[:,1].min()+local[:,1].max())/2,-local[:,2].min()]);placed=source.transform(np.c_[frame.T,shift]);part=write(output/(name+'.stl'),placed);tv=np.array(part['vertices']);tf=np.array(part['faces']);tri=tv[tf];mask=np.all(tri[:,:,2]<=.05,axis=1);bed_area=float(np.linalg.norm(np.cross(tri[mask,1]-tri[mask,0],tri[mask,2]-tri[mask,0]),axis=1).sum()/2)
            raw=(output/(name+'.stl')).read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);mesh=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(mesh.reshape(-1,3),axis=0,return_inverse=True);re=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)));relative=abs(re.volume()-source.volume())/max(1,abs(source.volume()));assert re.status()==mf.Error.NoError and relative<1e-4
            row={'muscle':fid,'part':name,'source_OBJ_sha256':hashlib.sha256((folder/(name+'.obj')).read_bytes()).hexdigest(),'print_STL_sha256':hashlib.sha256(raw).hexdigest(),'world_to_print_rotation_rows':frame.T.tolist(),'world_to_print_translation_mm':shift.tolist(),'print_bbox_mm':list(placed.bounding_box()),'triangular_area_within0.05mm_bed_mm2':bed_area,'STL_relative_volume_error':relative,'slicing_verified':False,'support_verified':False,'physical_fit_verified':False};results.append(row)
        (output/'PRINT_README.md').write_text('# M5刚性端帽尺寸样件\n\n此目录只含端帽底座与球窝保持盖，不含TPU肌腹、金属球头/压接坯和M2紧固件。\n模型为独立打印坐标，不能与原装配坐标直接叠加。\n未验证切片、支撑、材料、公差或强度；端帽底座的球窝/耳部可能低于原端板平面，必须检查首层与悬空结构。\n记录实际材料、喷嘴/层高、支撑、孔处理与尺寸，再检查手动装配、转向、压接滑移和损伤。\n金属球头坯的STL仅为形状包络，不作为TPU可承载打印件。\n',encoding='utf8')
    (DEST/'print_transforms.json').write_text(json.dumps({'scope':'four rigid housing print-placement candidates, not slicing or manufacturing approval','parts':results},indent=2),encoding='utf8');print('M5 rigid housing print candidates',len(results))
if __name__=='__main__':run()
