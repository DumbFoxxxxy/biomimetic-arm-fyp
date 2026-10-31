"""MuJoCo comparison of identical source anatomy and two generated structures."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import json
import numpy as np,mujoco
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parents[1]/'exports/anatomical_concept'

def run():
    fid='FJ1486';inv=json.loads((OUT/'inventory.json').read_text(encoding='utf8'));rec=next(r for r in inv['models'] if r['id']==fid);bb=np.array(rec['raw_bounds_mm'])-np.array(inv['translation_subtracted_mm']);centre=bb.mean(0)/1000;panels=[]
    for folder,title in [(OUT,'Source anatomy'),(OUT/'lattice_screen','Clipped BCC grid'),(OUT/'surface_cage_screen','BCC + surface cage')]:
        m=_load_project_model(str(folder/'reference.xml'));d=mujoco.MjData(m);mujoco.mj_forward(m,d)
        for g in range(m.ngeom):
            if mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g)!=fid:m.geom_rgba[g,3]=0
            else:m.geom_rgba[g,:]=[.8,.12,.10,1]
        cam=mujoco.MjvCamera();cam.lookat[:]=centre;cam.distance=.39;cam.azimuth=-90;cam.elevation=0;renderer=mujoco.Renderer(m,800,500);renderer.update_scene(d,camera=cam);im=Image.fromarray(renderer.render());draw=ImageDraw.Draw(im);draw.rectangle([0,0,500,52],fill='white');draw.text((15,15),title,fill='black');draw.rectangle([0,745,500,800],fill='white');draw.text((15,755),'FJ1486 | source shape retained',fill='black');draw.text((15,775),'No fiber/wire/actuation validation',fill='black');panels.append(im);renderer.close()
    canvas=Image.new('RGB',(1500,800),'white')
    for i,im in enumerate(panels):canvas.paste(im,(i*500,0))
    canvas.save(OUT/'lattice_detail_comparison.png');print(OUT/'lattice_detail_comparison.png')

if __name__=='__main__':run()
