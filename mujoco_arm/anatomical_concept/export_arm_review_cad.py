"""Combine the selected SC/AC/GH/E solids without duplicated superseded bones."""
from export_SC_cad import export_cad,ROOT
import json,hashlib
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection'
DEST=BASE/'AC_compact/GH_connection/elbow_connection'

def read(p):return json.loads(p.read_text(encoding='utf8'))

def run():
    selected={};layers=[(BASE,{'clavicle_SC_ball':'FJ3362','manubrium_socket_base':'FJ3290'}),(BASE/'AC_compact',{'clavicle_SC_AC_balls':'FJ3362','scapula_AC_socket':'FJ3384'}),(BASE/'AC_compact/GH_connection',None),(DEST,None)]
    for folder,aliases in layers:
        assembly=read(folder/'assembly.json');aliases=aliases if aliases is not None else assembly['bone_aliases']
        cad=read(folder/'cad_roundtrip.json')['parts'];inputs={p['name']:p for p in read(folder/'cad_input.json')['parts']}
        for name in assembly['components']:
            check=next((p for p in cad if p['part']==name),None)
            if not check or not check['pass']:raise RuntimeError('Unverified STEP '+name)
            volume=inputs[name]['volume_mm3']
            if abs(check['input_volume_mm3']-volume)>1e-6*max(1,volume):raise RuntimeError('Stale CAD evidence '+name)
            path=folder/'cad'/(name+'.step')
            selected[aliases.get(name,name)]={'part':name,'source_STEP':str(path.relative_to(ROOT)),'source_volume_mm3':volume,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    builder=export_cad.BRep_Builder();compound=export_cad.TopoDS_Compound();builder.MakeCompound(compound)
    for name,row in selected.items():
        reader=export_cad.STEPControl_Reader()
        if reader.ReadFile(str(ROOT/row['source_STEP']))!=export_cad.IFSelect_RetDone:raise RuntimeError('Read failed '+name)
        reader.TransferRoots();builder.Add(compound,reader.OneShape())
    path=DEST/'cad/anatomical_arm_review.step';shape=export_cad.export(compound,path);iterator=export_cad.TopExp_Explorer(shape,export_cad.TopAbs_SOLID);count=0
    while iterator.More():count+=1;iterator.Next()
    expected=sum(r['source_volume_mm3'] for r in selected.values());relative=abs(export_cad.volume(shape)-expected)/max(1,expected);valid=export_cad.BRepCheck_Analyzer(shape).IsValid()
    result={'scope':'whole selected rigid anatomical skeleton candidate; no muscles/support mount or native CAD mates','format':'faceted BRep compound, separate solids; editable import, not native source feature history','units':'mm','solid_count':count,'expected':len(selected),'BRep_valid':valid,'relative_volume_error':relative,'pass':valid and count==len(selected)==42 and relative<1e-4,'frozen':False,'parts':selected,'STEP':str(path.relative_to(ROOT)),'STEP_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    (DEST/'whole_arm_STEP_check.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    if not result['pass']:raise RuntimeError('Whole arm STEP mismatch')
    print('whole arm',count,'solids valid',valid,'relative volume error',relative,flush=True)

if __name__=='__main__':run()
