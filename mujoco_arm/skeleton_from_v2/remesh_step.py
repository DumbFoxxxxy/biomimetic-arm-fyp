"""Re-export the STEP BRep faces to triangles for an independent comparison."""
from export_step import *
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopLoc import TopLoc_Location
from OCP.TopAbs import TopAbs_FACE,TopAbs_REVERSED
def run():
    data=json.loads((OUT/'cad_input.json').read_text());meshes={}
    for part in data['parts']:
        reader=STEPControl_Reader();status=reader.ReadFile(str(OUT/'cad'/(part['name']+'.step')))
        if status!=IFSelect_RetDone:raise RuntimeError('Cannot read STEP')
        reader.TransferRoots();shape=reader.OneShape();BRepMesh_IncrementalMesh(shape,.005)
        exp=TopExp_Explorer(shape,TopAbs_FACE);triangles=[]
        while exp.More():
            face=TopoDS.Face_s(exp.Current());loc=TopLoc_Location();mesh=BRep_Tool.Triangulation_s(face,loc)
            if mesh is None:raise RuntimeError('Unmeshed face')
            transform=loc.Transformation()
            for i in range(1,mesh.NbTriangles()+1):
                indices=list(mesh.Triangle(i).Get())
                if face.Orientation()==TopAbs_REVERSED:indices[1],indices[2]=indices[2],indices[1]
                pts=[]
                for j in indices:
                    p=mesh.Node(j).Transformed(transform);pts.append([p.X(),p.Y(),p.Z()])
                triangles.append(pts)
            exp.Next()
        meshes[part['name']]=triangles;print(part['name'],len(triangles),flush=True)
    (OUT/'cad_roundtrip_meshes.json').write_text(json.dumps({'input_cache_sha256':data['input_cache_sha256'],'meshes_design_mm':meshes}))
if __name__=='__main__':run()
