"""Current M4 neutral109-rigid compound, preserving exact source placements."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.BRep import BRep_Tool
import json,hashlib,numpy as np,sys
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'
def read(path):return json.loads(path.read_text(encoding='utf8'))
def run():
    variant=sys.argv[sys.argv.index('--variant')+1] if '--variant' in sys.argv else 'axis_negative_r15';folder=BASE/'side_lever_candidate'/variant
    assembly=read(folder/'assembly.json');parent=read(BASE/'mounted_arm_STEP_check.json');assert parent['pass']
    inputs={p['name']:p for p in read(folder/'cad_input.json')['parts']};checks={p['part']:p for p in read(folder/'cad_roundtrip.json')['parts']}
    builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);parts={}
    for r in assembly['instances']:
        if r['kind']=='TPU':continue
        key=r['key'];name='forearm_spotfaced' if key=='FJ3391' else key
        if name in inputs:
            check=checks[name];assert check['pass'] and abs(check['input_volume_mm3']-inputs[name]['volume_mm3'])<1e-5
            path=folder/'cad'/(name+'.step');expected_volume=inputs[name]['volume_mm3']
        else:
            record=parent['parts'][key];path=ROOT/record['source_STEP'];expected_volume=record['source_volume_mm3']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==record['STEP_sha256']
            assert hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest()==record['OBJ_sha256']
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape();assert cad.BRepCheck_Analyzer(shape).IsValid()
        vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
        while it.More():
            p=BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
        vertices=np.unique(np.array(vertices),axis=0)
        expected=np.array([[float(x) for x in line.split()[1:4]] for line in (ROOT/r['file']).read_text(encoding='utf8').splitlines() if line.startswith('v ')])+np.array(r['translation_display_mm'])
        def directed(a,b):return max(float(np.sqrt(((a[i:i+128,None]-b[None])**2).sum(2)).min(1).max()) for i in range(0,len(a),128))
        error=max(directed(vertices,expected),directed(expected,vertices));assert error<.001,(key,error)
        assert abs(cad.volume(shape)-expected_volume)<1e-4*max(1,expected_volume)
        builder.Add(compound,shape)
        parts[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_OBJ':r['file'],'OBJ_sha256':hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest(),'translation_display_mm':r['translation_display_mm'],'source_volume_mm3':expected_volume,'STEP_vs_OBJ_vertex_error_mm':error,'owner':r['owner']}
    assert len(parts)==109
    path=folder/'cad/side_lever_arm_review.step';shape=cad.export(compound,path);it=cad.TopExp_Explorer(shape,cad.TopAbs_SOLID);count=0
    while it.More():count+=1;it.Next()
    expected=sum(r['source_volume_mm3'] for r in parts.values());relative=abs(cad.volume(shape)-expected)/expected;valid=cad.BRepCheck_Analyzer(shape).IsValid();assert count==109 and valid and relative<1e-4
    result={'scope':'M4 neutral109 rigid pieces, including20 unconnected distal-end pieces; twoTPU excluded. No whole muscle/drive validation.','format':'independent faceted BRep solids in mm; no native source feature history or assembly mates','solid_count':count,'relative_volume_error':relative,'BRep_valid':valid,'assembly_sha256':hashlib.sha256((folder/'assembly.json').read_bytes()).hexdigest(),'STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'parts':parts,'pass':True,'frozen':False}
    (folder/'side_lever_arm_STEP_check.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M4 wholeCAD',count,'valid',valid,'volume error',relative,flush=True)
if __name__=='__main__':run()
