"""105-piece neutral CAD of direct-muscle interfaces, two TPU meshes separate."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.gp import gp_Trsf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
import json,numpy as np,hashlib,sys
OUT=ROOT/'mujoco_arm/exports/anatomical_concept';BASE=OUT/'elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15';SHIFTED='--shifted' in sys.argv;DEST=OUT/'direct_muscle_interfaces'/('mechanical_yoke_shifted' if SHIFTED else 'mechanical_yoke_zero')
if '--access-relief' in sys.argv:DEST=OUT/'direct_muscle_interfaces/yoke_access_relief_study'
if '--extensor-reanchor' in sys.argv:DEST=OUT/'direct_muscle_interfaces/mechanical_yoke_reanchored_pair'

def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    a=read(DEST/'assembly.json');parent=read(BASE/'side_lever_arm_STEP_check.json');native={r['name']:r for r in read(DEST/'cad_input.json')['parts']};checks={r['part']:r for r in read(DEST/'cad_roundtrip.json')['parts']};assert all(r['pass'] for r in checks.values());delta=np.array(read(OUT/'direct_muscle_interfaces/flexor_attachment_search.json')['selected']['delta_display_mm']) if SHIFTED else np.zeros(3);builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);records={}
    for r in a['instances']:
        if r['kind']=='TPU':continue
        key=r['key'];name='direct_forearm_seated' if key=='FJ3391' else key;t=np.zeros(3)
        if name in native:path=DEST/'cad'/(name+'.step');volume=native[name]['volume_mm3'];assert checks[name]['pass']
        else:
            p=parent['parts'][key];path=ROOT/p['source_STEP'];assert digest(path)==p['STEP_sha256'];volume=p['source_volume_mm3']
            if SHIFTED and (key.startswith('FJ1486_distal_') or key=='FJ1486_active_wire_crimp_barrel'):t=delta
            if '--extensor-reanchor' in sys.argv and (key.startswith('FJ1477_distal_') or key=='FJ1477_active_wire_crimp_barrel'):t=np.array(read(OUT/'direct_muscle_interfaces/extensor_reanchor_extended/endpoint_search.json')['selected']['delta_display_mm'])
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape()
        if np.linalg.norm(t)>0:tr=gp_Trsf();tr.SetValues(*np.c_[np.eye(3),t].reshape(-1));shape=BRepBuilderAPI_Transform(shape,tr,True).Shape()
        assert cad.BRepCheck_Analyzer(shape).IsValid();vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
        while it.More():p=cad.BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
        vertices=np.unique(vertices,axis=0);expected=np.array([[float(x) for x in line.split()[1:4]] for line in (ROOT/r['file']).read_text(encoding='utf8').splitlines() if line.startswith('v ')])+np.array(r['translation_display_mm'])
        def directed(x,y):return max(float(np.sqrt(((x[i:i+128,None]-y[None])**2).sum(2)).min(1).max()) for i in range(0,len(x),128))
        error=max(directed(vertices,expected),directed(expected,vertices));assert error<.001,(key,error);assert abs(cad.volume(shape)-volume)<.0001*max(1,volume);builder.Add(compound,shape);records[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_OBJ':r['file'],'OBJ_sha256':digest(ROOT/r['file']),'source_volume_mm3':volume,'additional_STEP_translation_mm':t.tolist(),'source_OBJ_translation_mm':r['translation_display_mm'],'STEP_vs_OBJ_vertex_error_mm':error,'owner':r['owner'],'nominal_packing_only':r.get('nominal_packing_only',False)}
    path=DEST/'cad/direct_muscle_arm_review.step';loaded=cad.export(compound,path);it=cad.TopExp_Explorer(loaded,cad.TopAbs_SOLID);count=0
    while it.More():count+=1;it.Next()
    relative=abs(cad.volume(loaded)-sum(r['source_volume_mm3'] for r in records.values()))/sum(r['source_volume_mm3'] for r in records.values());assert count==len(records)==105 and cad.BRepCheck_Analyzer(loaded).IsValid() and relative<.0001
    (DEST/'whole_arm_STEP_check.json').write_text(json.dumps({'scope':'105 neutral rigid solids, including2uncrimped nominal stop poses; two TPU meshes separate, no muscle drive/strength proof','format':'faceted editable BRep, not native feature history/mates','solid_count':count,'valid':True,'relative_volume_error':relative,'STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_assembly_sha256':digest(DEST/'assembly.json'),'parts':records,'frozen':False},indent=2),encoding='utf8');print('direct-arm105 CAD valid, max sourcevertex error',max(r['STEP_vs_OBJ_vertex_error_mm'] for r in records.values()),flush=True)

if __name__=='__main__':run()
