"""27-part M8 conditional neutral interface CAD with explicit joint poses."""
from export_SC_cad import export_cad as cad,ROOT
from OCP.gp import gp_Trsf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
import json,numpy as np,hashlib
DEST=ROOT/'mujoco_arm/exports/anatomical_concept/forearm_receivers';MOUNT=DEST/'rephased_anchors';BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate/side_lever_candidate/axis_negative_r15'
YAW_PARTS={'rotor_yoke','pitch_sleeve','pitch_M2x12','pitch_front_washer','pitch_back_washer','pitch_nut','pitch_crimp_blank'}

def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    probe=read(MOUNT/'mounted_direction_screen.json');assert probe['first_rigid_limit'] is None and not probe['neutral_contacts'];definition=read(BASE/'design.json');axis=np.array(definition['joint_axis']);down=np.array(definition['zero_forearm_direction']);F=np.stack([down,np.cross(-axis,down),-axis],axis=1);builder=cad.BRep_Builder();compound=cad.TopoDS_Compound();builder.MakeCompound(compound);records=[]
    requests=[('lever_receiving_spotfaces',MOUNT/'cad/lever_receiving_spotfaces.step',MOUNT/'lever_receiving_spotfaces.obj',np.eye(3),np.zeros(3))]
    for fid,label in [('FJ1486','flexor'),('FJ1477','extensor')]:
        d=read(DEST/fid/'design.json');metrics=probe['poses'][0]['modules'][fid];yaw=np.radians(metrics['yaw_relative_lever_deg']);pitch=np.radians(metrics['pitch_deg']);Rz=np.array([[np.cos(yaw),-np.sin(yaw),0],[np.sin(yaw),np.cos(yaw),0],[0,0,1]]);Rp=np.array([[np.cos(pitch),0,-np.sin(pitch)],[0,1,0],[np.sin(pitch),0,np.cos(pitch)]]);H=np.array([8,0,5.5+d['standoff_mm']]);P=np.array(probe['candidate_eyelets_display_mm'][label])
        for key in d['parts']:
            R=np.eye(3);t=np.zeros(3)
            if key=='pitch_crimp_blank':R=Rp;t=H-Rp@H
            if key in YAW_PARTS:R=Rz@R;t=Rz@t
            requests.append((fid+'_'+key,DEST/fid/'cad'/(key+'.step'),MOUNT/(fid+'_'+key+'_neutral.obj'),F@R,P+F@t))
    for name,path,obj,R,t in requests:
        reader=cad.STEPControl_Reader();assert reader.ReadFile(str(path))==cad.IFSelect_RetDone;reader.TransferRoots();tr=gp_Trsf();tr.SetValues(*np.c_[R,t].reshape(-1));shape=BRepBuilderAPI_Transform(reader.OneShape(),tr,True).Shape();assert cad.BRepCheck_Analyzer(shape).IsValid();vertices=[];it=cad.TopExp_Explorer(shape,cad.TopAbs_VERTEX)
        while it.More():p=cad.BRep_Tool.Pnt_s(cad.TopoDS.Vertex_s(it.Current()));vertices.append([p.X(),p.Y(),p.Z()]);it.Next()
        vertices=np.unique(vertices,axis=0);expected=np.array([[float(x) for x in line.split()[1:4]] for line in obj.read_text(encoding='utf8').splitlines() if line.startswith('v ')])
        def directed(a,b):return max(float(np.sqrt(((a[i:i+128,None]-b[None])**2).sum(2)).min(1).max()) for i in range(0,len(a),128))
        error=max(directed(vertices,expected),directed(expected,vertices));assert error<.001,(name,error);builder.Add(compound,shape);records.append({'part':name,'source_STEP':str(path.relative_to(ROOT)),'STEP_sha256':digest(path),'source_OBJ':str(obj.relative_to(ROOT)),'OBJ_sha256':digest(obj),'rotation':R.tolist(),'translation_mm':t.tolist(),'source_vs_CAD_vertex_error_mm':error,'volume_mm3':cad.volume(shape),'valid':True})
    path=MOUNT/'cad/receiving_interfaces_review.step';loaded=cad.export(compound,path);it=cad.TopExp_Explorer(loaded,cad.TopAbs_SOLID);count=0
    while it.More():count+=1;it.Next()
    error=abs(cad.volume(loaded)-sum(r['volume_mm3'] for r in records))/sum(r['volume_mm3'] for r in records);assert count==27 and cad.BRepCheck_Analyzer(loaded).IsValid() and error<.0001
    (MOUNT/'assembly_STEP_check.json').write_text(json.dumps({'scope':'M8 conditional neutral 27-piece receiving interfaces only, poses reference sheave centres, not completed cable transmission or arm acceptance','solid_count':count,'valid':True,'relative_volume_error':error,'STEP_sha256':digest(path),'direction_screen_sha256':digest(MOUNT/'mounted_direction_screen.json'),'parts':records,'frozen':False},indent=2),encoding='utf8');print('M8 27-piece STEP valid, max vertex error',max(r['source_vs_CAD_vertex_error_mm'] for r in records),flush=True)

if __name__=='__main__':run()
