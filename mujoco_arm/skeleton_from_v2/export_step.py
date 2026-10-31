"""Create editable faceted BRep solids and verify a STEP round trip.

Run with Python3.10 and project-local OpenCascade dependencies. Faceted BRep
is explicit: this does not restore the original SolidWorks feature history.
"""
from pathlib import Path
import sys,json,math
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'mujoco_arm/exports/skeleton_from_v2'
sys.path.insert(0,str(ROOT/'tmp/skeleton_cad_vendor310'))
from OCP.gp import gp_Pnt
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakePolygon,BRepBuilderAPI_MakeFace,BRepBuilderAPI_Sewing,BRepBuilderAPI_MakeSolid
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SHELL,TopAbs_SOLID,TopAbs_VERTEX
from OCP.TopoDS import TopoDS,TopoDS_Compound
from OCP.BRep import BRep_Builder,BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepLib import BRepLib
from OCP.STEPControl import STEPControl_Writer,STEPControl_Reader,STEPControl_AsIs
from OCP.IFSelect import IFSelect_RetDone
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.ShapeFix import ShapeFix_Shape

def volume(shape):
    props=GProp_GProps();BRepGProp.VolumeProperties_s(shape,props);return props.Mass()
def bounds(shape):
    box=Bnd_Box();BRepBndLib.Add_s(shape,box);return list(box.Get())
def export(shape,path):
    writer=STEPControl_Writer()
    if writer.Transfer(shape,STEPControl_AsIs)!=IFSelect_RetDone:raise RuntimeError('STEP transfer failed')
    if writer.Write(str(path))!=IFSelect_RetDone:raise RuntimeError('STEP write failed')
    reader=STEPControl_Reader()
    if reader.ReadFile(str(path))!=IFSelect_RetDone:raise RuntimeError('STEP read failed')
    reader.TransferRoots();return reader.OneShape()
def run():
    data=json.loads((OUT/'cad_input.json').read_text());cad=OUT/'cad';cad.mkdir(exist_ok=True);rows=[];vertices={};builder=BRep_Builder();assembly=TopoDS_Compound();builder.MakeCompound(assembly)
    for part in data['parts']:
        points=[gp_Pnt(*p) for p in part['vertices_design_mm']];sew=BRepBuilderAPI_Sewing(1e-6)
        for face in part['faces']:
            polygon=BRepBuilderAPI_MakePolygon()
            for i in face:polygon.Add(points[i])
            polygon.Close();maker=BRepBuilderAPI_MakeFace(polygon.Wire())
            if not maker.IsDone():raise RuntimeError('Invalid triangle in '+part['name'])
            sew.Add(maker.Face())
        sew.Perform();shape=sew.SewedShape();solidmaker=BRepBuilderAPI_MakeSolid();explore=TopExp_Explorer(shape,TopAbs_SHELL);shells=0
        while explore.More():solidmaker.Add(TopoDS.Shell_s(explore.Current()));shells+=1;explore.Next()
        if shells<1:raise RuntimeError('No closed shells')
        solid=solidmaker.Solid();BRepLib.OrientClosedSolid_s(solid)
        healed=False
        if not BRepCheck_Analyzer(solid).IsValid():
            fixer=ShapeFix_Shape(solid);fixer.SetPrecision(1e-6);fixer.SetMaxTolerance(.001);fixer.Perform();solid=fixer.Shape();healed=True
        if not BRepCheck_Analyzer(solid).IsValid():raise RuntimeError('Invalid CAD solid after healing: '+part['name'])
        rebuilt=export(solid,cad/(part['name']+'.step'));valid=BRepCheck_Analyzer(rebuilt).IsValid();v=volume(rebuilt);box=bounds(rebuilt);originalbox=bounds(solid)
        relative=abs(v-part['volume_mm3'])/max(1,part['volume_mm3']);bboxerror=max(abs(a-b) for a,b in zip(box,originalbox));exp=TopExp_Explorer(rebuilt,TopAbs_SOLID);count=0
        while exp.More():count+=1;exp.Next()
        row={'bone':part['name'],'topology_healing_applied':healed,'input_volume_mm3':part['volume_mm3'],'roundtrip_volume_mm3':v,'relative_volume_error':relative,'bbox_roundtrip_error_mm':bboxerror,'valid_brep':valid,'solid_count':count,'closed_shell_count':shells,'pass':valid and count==1 and relative<1e-5 and bboxerror<1e-4}
        if not row['pass']:raise RuntimeError(str(row))
        pts=[];explorer=TopExp_Explorer(rebuilt,TopAbs_VERTEX)
        while explorer.More():
            p=BRep_Tool.Pnt_s(TopoDS.Vertex_s(explorer.Current()));pts.append([p.X(),p.Y(),p.Z()]);explorer.Next()
        vertices[part['name']]=pts
        rows.append(row);builder.Add(assembly,rebuilt);print(json.dumps(row),flush=True)
    combined=export(assembly,cad/'skeleton_zero_assembly.step')
    report={'input_cache_sha256':data['input_cache_sha256'],'format':'editable faceted BRep STEP, not original SolidWorks parametric feature history','coordinate_system':data['coordinate_system'],'parts':rows,'pass':all(r['pass'] for r in rows),'assembly_valid':BRepCheck_Analyzer(combined).IsValid()}
    (OUT/'cad_roundtrip.json').write_text(json.dumps(report,indent=2));print('CAD roundtrip complete',flush=True)
    (OUT/'cad_roundtrip_vertices.json').write_text(json.dumps(vertices))
if __name__=='__main__':run()
