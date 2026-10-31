"""Fuse verified parent forearm STEP with a separately valid rotor solid.

Avoid sewing Boolean-generated microscopic triangles in the final union mesh.
This remains a faceted BRep, not reconstructed anatomical feature history.
"""
from export_SC_cad import export_cad as cad,ROOT
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse,BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.gp import gp_Vec
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.TopoDS import TopoDS_Shell
from OCP.TopAbs import TopAbs_FACE
import numpy as np,json,hashlib
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/distal_rotor_space/pocket_candidate'
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
def read(path):return json.loads(path.read_text(encoding='utf8'))
def shape(path):
    reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone
    reader.TransferRoots();return reader.OneShape()
def run():
    if '--humerus-only' in __import__('sys').argv:
        import shutil
        folder=DEST/'humerus_only';folder.mkdir(exist_ok=True)
        part=next(p for p in read(DEST/'cad_input.json')['parts'] if p['name']=='humerus_rotor_pocket')
        (folder/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[part]}),encoding='utf8')
        cad.OUT=folder;cad.run();(DEST/'cad').mkdir(exist_ok=True)
        shutil.copy2(folder/'cad/humerus_rotor_pocket.step',DEST/'cad/humerus_rotor_pocket.step')
        shutil.copy2(folder/'cad_roundtrip.json',DEST/'cad_roundtrip.json');return
    # Recreate the same24-sided prismatic rings directly in BRep. Sewing the
    # union's microtriangles failed even before adding the anatomical forearm.
    selection=read(DEST/'selection.json')['geometry'];screen=read(DEST.parent/'space_screen.json')
    E=np.array(screen['centre_display_mm']);axis=np.array(screen['axis']);z=axis
    x=np.cross(z,[0,0,1] if abs(z[2])<.8 else [0,1,0]);x/=np.linalg.norm(x);y=np.cross(z,x)
    def cylinder(lo,hi,radius):
        polygon=cad.BRepBuilderAPI_MakePolygon()
        for angle in np.arange(24)*np.pi/12:
            point=E+axis*lo+radius*(x*np.cos(angle)+y*np.sin(angle));polygon.Add(cad.gp_Pnt(*point))
        polygon.Close();face=cad.BRepBuilderAPI_MakeFace(polygon.Wire()).Face()
        return BRepPrimAPI_MakePrism(face,gp_Vec(*(axis*(hi-lo)))).Shape()
    def fuse(a,b):
        op=BRepAlgoAPI_Fuse(a,b);op.Build();assert op.IsDone();return op.Shape()
    centre=selection['axis_station_mm'];radius=selection['nominal_steel_centre_radius_mm'];lo=centre-2;hi=centre+2
    start,end=sorted([np.sign(centre)*7.6,centre])
    rotor=fuse(fuse(cylinder(lo,lo+1.5,radius+1),cylinder(lo+1.5,hi-1.5,radius-.3)),cylinder(hi-1.5,hi,radius+1))
    rotor=fuse(rotor,cylinder(start,end,8))
    cut=BRepAlgoAPI_Cut(rotor,cylinder(min(lo,start)-1,max(hi,end)+1,6));cut.Build();assert cut.IsDone()
    rotor=cut.Shape();assert cad.BRepCheck_Analyzer(rotor).IsValid()
    parent=BASE/'cad/forearm_elbow_hub.step';addition=DEST/'rotor_addition/cad/rotor_addition.step'
    addition.parent.mkdir(exist_ok=True);rotor=cad.export(rotor,addition)
    rotor_input=read(DEST/'rotor_addition/cad_input.json')['parts'][0]
    rotor_error=abs(cad.volume(rotor)-rotor_input['volume_mm3'])
    print('native rotor volume',cad.volume(rotor),'mesh volume',rotor_input['volume_mm3'],'difference',rotor_error,flush=True)
    # The manufacturing union mesh was simplified at0.001mm; compare with the
    # same relative-volume threshold as other CAD roundtrips, not exact equality.
    assert rotor_error/max(1,abs(rotor_input['volume_mm3']))<1e-4
    operation=BRepAlgoAPI_Fuse(shape(parent),shape(addition));operation.Build()
    assert operation.IsDone();operation.SimplifyResult();result=operation.Shape()
    assert cad.BRepCheck_Analyzer(result).IsValid(),'native union invalid'
    path=DEST/'cad/forearm_with_rotor.step';reloaded=cad.export(result,path)
    assert cad.BRepCheck_Analyzer(reloaded).IsValid()
    if '--canonical-mesh' in __import__('sys').argv:
        from OCP.BRepMesh import BRepMesh_IncrementalMesh
        from OCP.StlAPI import StlAPI_Writer
        BRepMesh_IncrementalMesh(reloaded,.001,False,.05,False).Perform()
        writer=StlAPI_Writer();writer.ASCIIMode=False
        assert writer.Write(reloaded,str(DEST/'forearm_from_native.stl'))
        print('Canonical native-union mesh exported for independent topology/geometry check',flush=True)
        return
    source=next(p for p in read(DEST/'cad_input.json')['parts'] if p['name']=='forearm_with_rotor')
    vertices=np.array(source['vertices']);bounds=np.r_[vertices.min(0),vertices.max(0)]
    relative=abs(cad.volume(reloaded)-source['volume_mm3'])/max(1,abs(source['volume_mm3']))
    bbox_error=float(abs(np.array(cad.bounds(reloaded))-bounds).max())
    # Compare samples to the actual surface shell, not the filled solid, which
    # could give zero for vertices inside an incorrect oversized volume.
    builder=cad.BRep_Builder();shell=TopoDS_Shell();builder.MakeShell(shell)
    iterator=cad.TopExp_Explorer(reloaded,TopAbs_FACE)
    while iterator.More():builder.Add(shell,iterator.Current());iterator.Next()
    indices=np.unique(np.linspace(0,len(vertices)-1,min(512,len(vertices))).astype(int));maximum=0.;worst=None
    for i in indices:
        distance=BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(cad.gp_Pnt(*vertices[i])).Vertex(),shell)
        distance.Perform();assert distance.IsDone()
        if distance.Value()>maximum:
            maximum=distance.Value();nearest=distance.PointOnShape2(1)
            worst={'vertex_index':int(i),'point_mm':vertices[i].tolist(),'nearest_STEP_point_mm':[nearest.X(),nearest.Y(),nearest.Z()],'local_to_elbow_in_rotor_frame_mm':((vertices[i]-E)@np.stack([x,y,axis],1)).tolist()}
    row={'part':'forearm_with_rotor','construction':'native BRep Boolean union of separately verified source forearm and rotor; microscopic final union mesh was not sewn','valid':True,'input_volume_mm3':source['volume_mm3'],'input_vertex_count':len(vertices),'relative_volume_error':relative,'bbox_error_mm':bbox_error,'sampled_mesh_vertex_to_STEP_surface_max_mm':maximum,'surface_sample_count':len(indices),'surface_check_scope':'512 deterministic vertex samples; not exhaustive surface equivalence proof','source_STEP_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [parent,addition]},'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'pass':relative<1e-4 and bbox_error<.001 and maximum<.02}
    row['worst_surface_sample']=worst
    (DEST/'native_union_diagnostic.json').write_text(json.dumps(row,indent=2),encoding='utf8')
    assert row['pass'],row
    existing=read(DEST/'cad_roundtrip.json');hum=next(p for p in existing['parts'] if p['part']=='humerus_rotor_pocket');assert hum['pass']
    report={'scope':'two standalone candidate rigid CAD parts; no assembly, tendon or drive validation','parts':[hum,row],'cad_input_sha256':hashlib.sha256((DEST/'cad_input.json').read_bytes()).hexdigest(),'pass':True,'frozen':False}
    (DEST/'rotor_cad_roundtrip.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(row,indent=2),flush=True)
if __name__=='__main__':run()
