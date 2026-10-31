"""Combined review STEP of all SC candidate solids, preserving separate bodies."""
from export_SC_cad import *
import json

def run():
    out=export_cad.OUT;data=json.loads((out/'cad_input.json').read_text(encoding='utf8'));builder=export_cad.BRep_Builder();compound=export_cad.TopoDS_Compound();builder.MakeCompound(compound)
    for part in data['parts']:
        reader=export_cad.STEPControl_Reader()
        if reader.ReadFile(str(out/'cad'/(part['name']+'.step')))!=export_cad.IFSelect_RetDone:raise RuntimeError('STEP input failed')
        reader.TransferRoots();builder.Add(compound,reader.OneShape())
    label='elbow' if '--elbow' in sys.argv else 'GH' if '--GH' in sys.argv else 'AC' if '--AC' in sys.argv else 'SC'
    re=export_cad.export(compound,out/('cad/'+label+'_review_assembly.step'));explore=export_cad.TopExp_Explorer(re,export_cad.TopAbs_SOLID);count=0
    while explore.More():count+=1;explore.Next()
    rel=abs(export_cad.volume(re)-sum(p['volume_mm3'] for p in data['parts']))/max(1,sum(p['volume_mm3'] for p in data['parts']));valid=export_cad.BRepCheck_Analyzer(re).IsValid();report={'solid_count':count,'expected':len(data['parts']),'valid_BRep':valid,'relative_volume_error':rel,'pass':valid and count==len(data['parts']) and rel<1e-4,'format':'faceted BRep review compound, separate solids, not native assembly mates or source feature history'};(out/'assembly_STEP_check.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)

if __name__=='__main__':run()
