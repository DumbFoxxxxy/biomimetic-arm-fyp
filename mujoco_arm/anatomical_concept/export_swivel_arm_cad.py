"""M5 neutral129-rigid assembly plus11-piece individual terminal compounds."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.BRep import BRep_Tool
import json,hashlib,numpy as np,sys
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15'
DEST=BASE/('swivel_terminal_distal_load' if '--distal-load' in sys.argv else 'swivel_terminal_candidate')
def read(p):return json.loads(p.read_text(encoding='utf8'))
def run():
    assembly=read(DEST/'assembly.json');parent=read(BASE/'side_lever_arm_STEP_check.json');assert parent['pass']
    builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);per_module={};selected={};inputs={};checks={}
    for fid in ['FJ1486','FJ1477']:
        inputs[fid]={p['name']:p for p in read(DEST/fid/'cad_input.json')['parts']};checks[fid]={p['part']:p for p in read(DEST/fid/'cad_roundtrip.json')['parts']}
        per_module[fid]=cad.TopoDS_Compound();builder.MakeCompound(per_module[fid])
    for r in assembly['instances']:
        if r['kind']=='TPU':continue
        key=r['key'];fid=r.get('m5_module')
        if fid:
            name=r['m5_part'];p=inputs[fid][name];assert checks[fid][name]['pass'] and abs(checks[fid][name]['input_volume_mm3']-p['volume_mm3'])<1e-5
            path=DEST/fid/'cad'/(name+'.step');expected_volume=p['volume_mm3']
        else:
            row=parent['parts'][key];path=ROOT/row['source_STEP'];expected_volume=row['source_volume_mm3'];assert hashlib.sha256(path.read_bytes()).hexdigest()==row['STEP_sha256']
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape();assert cad.BRepCheck_Analyzer(shape).IsValid()
        vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
        while it.More():
            p=BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
        vertices=np.unique(np.array(vertices),axis=0);expected=np.array([[float(x) for x in line.split()[1:4]] for line in (ROOT/r['file']).read_text(encoding='utf8').splitlines() if line.startswith('v ')])+np.array(r['translation_display_mm'])
        def directed(a,b):return max(float(np.sqrt(((a[i:i+128,None]-b[None])**2).sum(2)).min(1).max()) for i in range(0,len(a),128))
        delta=max(directed(vertices,expected),directed(expected,vertices));assert delta<.001,(key,delta)
        assert abs(cad.volume(shape)-expected_volume)<1e-4*max(1,expected_volume)
        builder.Add(compound,shape)
        if fid:builder.Add(per_module[fid],shape)
        selected[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_OBJ':r['file'],'OBJ_sha256':hashlib.sha256((ROOT/r['file']).read_bytes()).hexdigest(),'translation_display_mm':r['translation_display_mm'],'source_volume_mm3':expected_volume,'STEP_vs_OBJ_vertex_error_mm':delta,'owner':r['owner'],'module':fid}
    def export_check(shape,path,expected_count,expected_volume):
        shape=cad.export(shape,path);it=cad.TopExp_Explorer(shape,cad.TopAbs_SOLID);count=0
        while it.More():count+=1;it.Next()
        relative=abs(cad.volume(shape)-expected_volume)/max(1,expected_volume);valid=cad.BRepCheck_Analyzer(shape).IsValid();assert count==expected_count and valid and relative<1e-4
        return {'solid_count':count,'BRep_valid':valid,'relative_volume_error':relative,'STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'pass':True}
    assert len(selected)==129
    (DEST/'cad').mkdir(exist_ok=True)
    whole=export_check(compound,DEST/'cad/swivel_arm_review.step',129,sum(r['source_volume_mm3'] for r in selected.values()));module_checks={}
    for fid,shape in per_module.items():
        module_checks[fid]=export_check(shape,DEST/fid/'cad/swivel_terminal_assembly.step',11,sum(r['source_volume_mm3'] for r in selected.values() if r['module']==fid))
    result={'scope':'M5 neutral packaging129 rigid pieces, not full installed compliant motion/strength/actuation','format':'separate faceted BRep solids, mm, no native feature history/mates','whole':whole,'modules':module_checks,'parts':selected,'source_assembly_sha256':hashlib.sha256((DEST/'assembly.json').read_bytes()).hexdigest(),'frozen':False}
    (DEST/'whole_arm_STEP_check.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M5 neutralCAD129 plus11/11 endpoints valid, volume error',whole['relative_volume_error'],flush=True)
if __name__=='__main__':run()
