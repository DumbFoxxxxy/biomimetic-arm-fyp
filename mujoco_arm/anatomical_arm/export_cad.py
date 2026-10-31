"""STEP geometry round-trip with project-local Python 3.10 / OpenCascade.
Full original-resolution bones and lattice solids. May take several minutes.
"""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'tmp/skeleton_cad_vendor310'),str(ROOT/'mujoco_arm/skeleton_from_v2')]
from export_step import *
OUT=ROOT/'mujoco_arm/exports/anatomical_arm/A1'

def run():
    global OUT
    if '--elbow-only' in sys.argv:OUT=OUT/'elbow_connection'
    if '--mounts-only' in sys.argv:OUT=OUT/'elbow_connection/muscle_mounts'
    if '--fittings-only' in sys.argv:OUT=OUT/'elbow_connection/muscle_mounts/end_fittings'
    cad=OUT/'cad';cad.mkdir(exist_ok=True);data=json.loads((OUT/'cad_input.json').read_text(encoding='utf8'));rows=[]
    prior_path=OUT/'cad_roundtrip.json';prior=json.loads(prior_path.read_text(encoding='utf8'))['parts'] if prior_path.exists() else []
    for part in data['parts']:
        if '--changed-only' in sys.argv and part['name'] in ['clavicle','scapula','humerus','ulna','radius']:
            cached=next((r for r in prior if r['part']==part['name'] and r['pass']),None)
            if cached and (cad/(part['name']+'.step')).exists():rows.append(cached);print('unchanged original bone STEP retained',part['name'],flush=True);continue
        print('CAD start',part['name'],len(part['faces']),flush=True);sew=BRepBuilderAPI_Sewing(1e-7);points=[gp_Pnt(*p) for p in part['vertices']]
        for face in part['faces']:
            polygon=BRepBuilderAPI_MakePolygon()
            for i in face:polygon.Add(points[i])
            polygon.Close();maker=BRepBuilderAPI_MakeFace(polygon.Wire())
            if maker.IsDone():sew.Add(maker.Face())
        sew.Perform();shape=sew.SewedShape();maker=BRepBuilderAPI_MakeSolid();explore=TopExp_Explorer(shape,TopAbs_SHELL);count=0
        while explore.More():maker.Add(TopoDS.Shell_s(explore.Current()));count+=1;explore.Next()
        if not count:raise RuntimeError('No closed shell: '+part['name'])
        shape=maker.Solid();BRepLib.OrientClosedSolid_s(shape)
        if not BRepCheck_Analyzer(shape).IsValid():
            fixer=ShapeFix_Shape(shape);fixer.SetPrecision(1e-5);fixer.SetMaxTolerance(.001);fixer.Perform();shape=fixer.Shape()
        valid=BRepCheck_Analyzer(shape).IsValid()
        if not valid:raise RuntimeError('Invalid BRep '+part['name'])
        reloaded=export(shape,cad/(part['name']+'.step'));v=volume(reloaded);relative=abs(v-part['volume_mm3'])/max(1,abs(part['volume_mm3']));bb=max(abs(a-b) for a,b in zip(bounds(shape),bounds(reloaded)))
        row={'part':part['name'],'shell_count':count,'input_volume_mm3':part['volume_mm3'],'input_vertex_count':len(part['vertices']),'valid':BRepCheck_Analyzer(reloaded).IsValid(),'relative_volume_error':relative,'bbox_error_mm':bb,'pass':relative<1e-4 and bb<1e-3 and BRepCheck_Analyzer(reloaded).IsValid()};rows.append(row)
        (OUT/'cad_roundtrip.json').write_text(json.dumps({'parts':rows,'format':'editable faceted BRep, not reconstructed native feature tree'},indent=2),encoding='utf8');print(row,flush=True)
        if not row['pass']:raise RuntimeError('STEP mismatch')

if __name__=='__main__':run()
