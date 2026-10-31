"""115-solid nominal neutral STEP, source-to-CAD geometry checked afresh."""
from export_direct_yoke_arm import cad,ROOT,read,digest,gp_Trsf,BRepBuilderAPI_Transform
import numpy as np,json,sys
PARENT=ROOT/'mujoco_arm/exports/anatomical_concept/direct_muscle_interfaces/mechanical_yoke_reanchored_pair'
DEST=PARENT/'paired_swivel_support'
EXPECTED_COUNT,NEW_COUNT=115,11
if '--keepers' in sys.argv:
    PARENT=DEST;DEST=PARENT/'captive_keepers';EXPECTED_COUNT,NEW_COUNT=125,15

def run():
    a=read(DEST/'assembly.json');parent=read(PARENT/'whole_arm_STEP_check.json');native={r['name']:r for r in read(DEST/'cad_input.json')['parts']};checks={r['part']:r for r in read(DEST/'cad_roundtrip.json')['parts']};assert len(native)==len(checks)==NEW_COUNT and all(r['pass'] for r in checks.values());builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);records={};start=digest(DEST/'assembly.json')
    for r in a['instances']:
        if r['kind']=='TPU':continue
        key=r['key'];translation=np.zeros(3)
        if key in native:path=DEST/'cad'/(key+'.step');volume=native[key]['volume_mm3']
        else:
            p=parent['parts'][key];path=ROOT/p['source_STEP'];assert digest(path)==p['STEP_sha256'];translation=np.array(p['additional_STEP_translation_mm']);volume=p['source_volume_mm3']
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape()
        if np.linalg.norm(translation)>0:
            tr=gp_Trsf();tr.SetValues(*np.c_[np.eye(3),translation].reshape(-1));shape=BRepBuilderAPI_Transform(shape,tr,True).Shape()
        assert cad.BRepCheck_Analyzer(shape).IsValid();vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
        while it.More():p=cad.BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
        vertices=np.unique(vertices,axis=0);expected=np.array([[float(x) for x in line.split()[1:4]] for line in (ROOT/r['file']).read_text(encoding='utf8').splitlines() if line.startswith('v ')])+np.array(r['translation_display_mm'])
        def directed(x,y):return max(float(np.sqrt(((x[i:i+128,None]-y[None])**2).sum(2)).min(1).max()) for i in range(0,len(x),128))
        error=max(directed(vertices,expected),directed(expected,vertices));assert error<.001,(key,error);assert abs(cad.volume(shape)-volume)<1e-4*max(1,volume);builder.Add(compound,shape);records[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_OBJ':r['file'],'OBJ_sha256':digest(ROOT/r['file']),'source_volume_mm3':volume,'additional_STEP_translation_mm':translation.tolist(),'source_OBJ_translation_mm':r['translation_display_mm'],'STEP_vs_OBJ_vertex_error_mm':error,'owner':r['owner'],'nominal_packing_only':r.get('nominal_packing_only',False)}
    path=DEST/'cad/paired_swivel_arm_review.step';loaded=cad.export(compound,path);it=cad.TopExp_Explorer(loaded,cad.TopAbs_SOLID);count=0
    while it.More():count+=1;it.Next()
    relative=abs(cad.volume(loaded)-sum(r['source_volume_mm3'] for r in records.values()))/sum(r['source_volume_mm3'] for r in records.values());assert count==len(records)==EXPECTED_COUNT and cad.BRepCheck_Analyzer(loaded).IsValid() and relative<1e-4 and start==digest(DEST/'assembly.json')
    (DEST/'whole_arm_STEP_check.json').write_text(json.dumps({'scope':f'{count} nominal neutral rigid solids,including2uncrimped nominal stops; TPU meshes separate. All cap swivel angles held0inCAD, conditional path usesdifferent prescribedposes; no equilibrium/drive proof','format':'faceted BRep, not native feature history orassemblymates','solid_count':count,'valid':True,'relative_volume_error':relative,'STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_assembly_sha256':start,'parts':records,'frozen':False},indent=2),encoding='utf8');print('paired-arm',count,'CADvalid maxsourcevertexerror',max(r['STEP_vs_OBJ_vertex_error_mm'] for r in records.values()),flush=True)

if __name__=='__main__':run()
