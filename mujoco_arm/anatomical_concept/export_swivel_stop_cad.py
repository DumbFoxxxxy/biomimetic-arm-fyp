"""Exact native Boolean shoulders avoid invalid sewn micro-triangle seams."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse,BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCP.gp import gp_Vec
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.StlAPI import StlAPI_Writer
import json,numpy as np,hashlib,shutil,subprocess
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15/swivel_terminal_distal_load'
def read(path):return json.loads(path.read_text(encoding='utf8'))
def run():
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=read(folder/'design.json');data=read(folder/'cad_input.json');other=folder/'other_hardware_source';other.mkdir(exist_ok=True);(other/'cad_input.json').write_text(json.dumps({'units':'mm','parts':[p for p in data['parts'] if p['name']!='distal_socket_base']}),encoding='utf8');cad.OUT=other;cad.run()
        cad.OUT=folder/'native_base_source';cad.run();path=cad.OUT/'cad/distal_socket_base.step';reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();shape=reader.OneShape();source_sha=hashlib.sha256(path.read_bytes()).hexdigest();u=np.array(design['selected']['neutral_terminal_axis']);ref=[0,0,1] if abs(u[2])<.8 else [0,1,0];x=np.cross(u,ref);x/=np.linalg.norm(x);y=np.cross(u,x)
        def cylinder(p,lo,hi,radius):
            poly=cad.BRepBuilderAPI_MakePolygon()
            for a in np.arange(24)*np.pi/12:poly.Add(cad.gp_Pnt(*(p+u*lo+radius*(x*np.cos(a)+y*np.sin(a)))))
            poly.Close();face=cad.BRepBuilderAPI_MakeFace(poly.Wire()).Face();return BRepPrimAPI_MakePrism(face,gp_Vec(*(u*(hi-lo)))).Shape()
        for centre in design['selected']['M2_centres_display_mm']:
            p=np.array(centre);cut=BRepAlgoAPI_Cut(cylinder(p,-.15,.15,2.5),cylinder(p,-1,1,1.1));cut.Build();assert cut.IsDone();fuse=BRepAlgoAPI_Fuse(shape,cut.Shape());fuse.SetFuzzyValue(1e-6);fuse.Build();assert fuse.IsDone();shape=fuse.Shape()
        assert cad.BRepCheck_Analyzer(shape).IsValid();(folder/'cad').mkdir(exist_ok=True);shape=cad.export(shape,folder/'cad/distal_socket_base.step');BRepMesh_IncrementalMesh(shape,.001,False,.05,False).Perform();writer=StlAPI_Writer();writer.ASCIIMode=False;assert writer.Write(shape,str(folder/'canonical_stops.stl'))
        subprocess.run(['D:/MujocoEnv/Scripts/python.exe',str(ROOT/'mujoco_arm/anatomical_concept/import_swivel_stop_mesh.py'),'--distal-load','--muscle',fid],cwd=ROOT,check=True)
        source=next(p for p in read(folder/'cad_input.json')['parts'] if p['name']=='distal_socket_base');relative=abs(cad.volume(shape)-source['volume_mm3'])/source['volume_mm3'];v=np.array(source['vertices']);bbox=np.r_[v.min(0),v.max(0)];bbox_error=float(abs(np.array(cad.bounds(shape))-bbox).max());assert relative<1e-4 and bbox_error<.001
        row={'part':'distal_socket_base','shell_count':1,'input_volume_mm3':source['volume_mm3'],'input_vertex_count':len(v),'valid':True,'relative_volume_error':relative,'bbox_error_mm':bbox_error,'construction':'native Boolean shoulders on verified cap, then CAD-derived canonical mesh','unstopped_source_STEP_sha256':source_sha,'pass':True}
        checks=read(other/'cad_roundtrip.json')['parts']
        for p in checks:shutil.copy2(other/'cad'/(p['part']+'.step'),folder/'cad'/(p['part']+'.step'))
        (folder/'cad_roundtrip.json').write_text(json.dumps({'parts':[row,*checks],'format':'faceted BRep with native Boolean stop shoulders; no native parametric history'},indent=2),encoding='utf8');print(fid,'native stops CAD valid',relative,bbox_error,flush=True)
if __name__=='__main__':run()
