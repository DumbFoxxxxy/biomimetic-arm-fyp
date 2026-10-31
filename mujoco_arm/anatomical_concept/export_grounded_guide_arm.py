"""M7 full neutral rigid compound, with explicit unchanged M5 proof reuse."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.gp import gp_Trsf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
import json,numpy as np,hashlib,sys
OUT=ROOT/'mujoco_arm/exports/anatomical_concept';GUIDES=OUT/'elbow_guides'
DEST=GUIDES/('ground_frame_swept' if '--swept-root' in sys.argv else 'ground_frame_low_root' if '--low-root' in sys.argv else 'ground_frame_candidate')
PARENT=OUT/'elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15/swivel_terminal_distal_load'

def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    a=read(DEST/'assembly.json');parent=read(PARENT/'whole_arm_STEP_check.json');native={p['name']:p for p in read(DEST/'cad_input.json')['parts']};checks={p['part']:p for p in read(DEST/'cad_roundtrip.json')['parts']};template={p['name']:p for p in read(GUIDES/'template/cad_input.json')['parts']};selection=read(GUIDES/'paired_placement/selection.json')['placements'];builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);attachments=cad.TopoDS_Compound();builder.MakeCompound(attachments);records={}
    for r in a['instances']:
        if r['kind']=='TPU':continue
        key=r['key'];new=True;transform=None;expected_obj=ROOT/r['file']
        if key in native:
            assert checks[key]['pass'];path=DEST/'cad'/(key+'.step');expected_volume=native[key]['volume_mm3']
        elif r.get('grounded_guide_candidate'):
            fid=key[:6];name=key[7:];assert name in template;path=GUIDES/'template/cad'/(name+'.step');expected_volume=template[name]['volume_mm3'];s=selection[fid];F=np.array(s['frame_columns']);Q=np.array(s['centre_display_mm']);transform=gp_Trsf();transform.SetValues(*(np.c_[F,Q].reshape(-1)))
        else:
            new=False;p=parent['parts'][key];path=ROOT/p['source_STEP'];expected_volume=p['source_volume_mm3'];assert digest(path)==p['STEP_sha256'];assert digest(expected_obj)==p['OBJ_sha256'];assert r['translation_display_mm']==p['translation_display_mm']
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape()
        if transform is not None:shape=BRepBuilderAPI_Transform(shape,transform,True).Shape()
        assert cad.BRepCheck_Analyzer(shape).IsValid();delta=None
        if new:
            vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
            while it.More():p=cad.BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
            vertices=np.unique(vertices,axis=0);expected=np.array([[float(x) for x in line.split()[1:4]] for line in expected_obj.read_text(encoding='utf8').splitlines() if line.startswith('v ')])+np.array(r['translation_display_mm'])
            def directed(x,y):return max(float(np.sqrt(((x[i:i+128,None]-y[None])**2).sum(2)).min(1).max()) for i in range(0,len(x),128))
            delta=max(directed(vertices,expected),directed(expected,vertices));assert delta<.001,(key,delta)
            builder.Add(attachments,shape)
        assert abs(cad.volume(shape)-expected_volume)<1e-4*max(1,expected_volume)
        builder.Add(compound,shape);records[key]={'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_OBJ':r['file'],'OBJ_sha256':digest(expected_obj),'source_volume_mm3':expected_volume,'vertex_error_mm_fresh':delta,'unchanged_parent_vertex_proof_reused':not new,'owner':r['owner']}
    assert len(records)==151
    def export(shape,filename,count,volume):
        path=DEST/'cad'/filename;loaded=cad.export(shape,path);it=cad.TopExp_Explorer(loaded,cad.TopAbs_SOLID);actual=0
        while it.More():actual+=1;it.Next()
        error=abs(cad.volume(loaded)-volume)/max(1,volume);valid=cad.BRepCheck_Analyzer(loaded).IsValid();assert actual==count and valid and error<1e-4
        return {'path':str(path.relative_to(ROOT)),'solid_count':actual,'valid':valid,'relative_volume_error':error,'sha256':digest(path),'pass':True}
    whole=export(compound,'grounded_guides_arm_review.step',151,sum(r['source_volume_mm3'] for r in records.values()));changed=[r for r in records.values() if not r['unchanged_parent_vertex_proof_reused']];sub=export(attachments,'guide_frame_subassemblies.step',24,sum(r['source_volume_mm3'] for r in changed));result={'scope':'151 rigid neutral packing solids, including free distal/crimp blanks; no compliant motion or force proof','format':'faceted BRep solids, no feature history or assembly mates','whole':whole,'changed_subassemblies':sub,'parts':records,'fresh_changed_vertex_checks':len(changed),'unchanged_parent_proofs_reused':127,'parent_manifest_sha256':digest(PARENT/'whole_arm_STEP_check.json'),'assembly_sha256':digest(DEST/'assembly.json'),'frozen':False};(DEST/'whole_arm_STEP_check.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('M7 151 rigid CAD,24 fresh source checks,127 unchanged proofs reused',flush=True)

if __name__=='__main__':run()
