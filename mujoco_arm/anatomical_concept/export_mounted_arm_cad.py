"""98 rigid neutral manufactured instances; no soft body or motion claim."""
from export_SC_cad import export_cad,ROOT
import json,hashlib
import numpy as np

DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'
SKELETON=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
def read(path):return json.loads(path.read_text(encoding='utf8'))
def run():
    module=DEST/'rigid_module_parts'
    if '--prepare-only-cad' in __import__('sys').argv:
        export_cad.OUT=module;export_cad.run();return
    assembly=read(DEST/'assembly.json');old=read(SKELETON/'whole_arm_STEP_check.json')['parts']
    selected={};inputs={};checks={}
    for folder in [DEST,module]:
        inputs[folder]={p['name']:p for p in read(folder/'cad_input.json')['parts']}
        checks[folder]={p['part']:p for p in read(folder/'cad_roundtrip.json')['parts']}
    builder=export_cad.BRep_Builder();compound=export_cad.TopoDS_Compound();builder.MakeCompound(compound)
    for r in assembly['instances']:
        if r['kind']=='TPU':continue
        key=r['key']
        if r.get('generated') or r.get('source_module'):
            folder=DEST if r.get('generated') else module
            name=__import__('pathlib').Path(r['file']).stem if r.get('generated') else key
            p=inputs[folder][name];check=checks[folder][name]
            assert check['pass'] and abs(check['input_volume_mm3']-p['volume_mm3'])<1e-5
            path=folder/'cad'/(name+'.step');expected_volume=p['volume_mm3']
        else:
            record=old[key];path=ROOT/record['source_STEP'];expected_volume=record['source_volume_mm3']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==record['source_sha256']
        reader=export_cad.STEPControl_Reader();assert reader.ReadFile(str(path))==export_cad.IFSelect_RetDone
        reader.TransferRoots();shape=reader.OneShape();assert export_cad.BRepCheck_Analyzer(shape).IsValid()
        # Compare every STEP vertex in world coordinates to the exact assembly OBJ.
        obj=(ROOT/r['file']).read_text(encoding='utf8').splitlines()
        points=np.array([[float(v) for v in line.split()[1:4]] for line in obj if line.startswith('v ')])+np.array(r['translation_display_mm'])
        actual=[];iterator=export_cad.TopExp_Explorer(shape,export_cad.TopAbs_VERTEX)
        from OCP.BRep import BRep_Tool
        while iterator.More():
            pnt=BRep_Tool.Pnt_s(export_cad.TopoDS.Vertex_s(iterator.Current()));actual.append([pnt.X(),pnt.Y(),pnt.Z()]);iterator.Next()
        actual=np.unique(np.array(actual),axis=0)
        def directed(a,b):
            return max(float(np.sqrt(((a[i:i+128,None,:]-b[None,:,:])**2).sum(2)).min(1).max()) for i in range(0,len(a),128))
        delta=max(directed(actual,points),directed(points,actual));assert delta<.001,(key,delta)
        volume=export_cad.volume(shape);assert abs(volume-expected_volume)<1e-4*max(1,abs(expected_volume))
        builder.Add(compound,shape)
        selected[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_OBJ':r['file'],'OBJ_sha256':hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest(),'translation_display_mm':r['translation_display_mm'],'source_volume_mm3':expected_volume,'STEP_vs_OBJ_vertex_error_mm':delta,'owner':r['owner'],'distal_unvalidated':r['owner']=='distal_unvalidated'}
        print('verified',key,'coordinate error mm',delta,flush=True)
    assert len(selected)==98
    path=DEST/'cad/mounted_arm_rigid_review.step';shape=export_cad.export(compound,path)
    iterator=export_cad.TopExp_Explorer(shape,export_cad.TopAbs_SOLID);count=0
    while iterator.More():count+=1;iterator.Next()
    expected=sum(r['source_volume_mm3'] for r in selected.values());relative=abs(export_cad.volume(shape)-expected)/max(1,expected)
    valid=export_cad.BRepCheck_Analyzer(shape).IsValid();assert valid and count==98 and relative<1e-4
    result={'scope':'98 manufactured rigid neutral pieces including free distal ends; TPU2 excluded. No whole-module motion or drive claim.','format':'separate faceted BRep solids in mm, no native mates/history','solid_count':count,'BRep_valid':valid,'relative_volume_error':relative,'parts':selected,'assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest(),'STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'pass':True,'frozen':False}
    (DEST/'mounted_arm_STEP_check.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('Mounted neutral STEP',count,'solids, relative volume error',relative,flush=True)

if __name__=='__main__':run()
