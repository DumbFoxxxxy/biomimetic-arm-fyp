"""Orthographic technical views of exported triangle arrays (millimetres)."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys
import numpy as np
from PIL import Image, ImageDraw

OUT=Path(__file__).resolve().parents[1]/'exports/skeleton_from_v2'
COLORS=[(160,174,186),(191,177,130),(89,167,135),(219,172,73),(105,145,204)]
def render(src,dest):
    meshes=np.load(src);im=Image.new('RGB',(1500,900),'#f3f5f7');draw=ImageDraw.Draw(im)
    for panel,(title,az) in enumerate([('Front oblique',25),('Rear oblique',155),('Joint detail',-45)]):
        a=np.radians(az);right=np.array([np.cos(a),np.sin(a),0]);up=np.array([0,0,1]);depth=np.cross(right,up)
        center=np.array([80.,0.,185. if panel==2 else 0.]);scale=3.3 if panel==2 else 1.25
        polygons=[]
        for idx,name in enumerate(meshes.files):
            if panel==2 and name=='thorax_geom':continue
            for tri in meshes[name]:
                xy=np.column_stack([(tri-center)@right,-(tri-center)@up])*scale+[panel*500+250,450]
                normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=max(np.linalg.norm(normal),1e-12)
                light=.65+.35*abs(normal@np.array([.3,-.5,.8]));color=tuple(int(c*light) for c in COLORS[idx])
                polygons.append((tri.mean(0)@depth,[tuple(p) for p in xy],color))
        for _,poly,color in sorted(polygons,key=lambda p:p[0]):draw.polygon(poly,fill=color)
        draw.text((panel*500+20,20),title,fill='black')
    im.save(dest)
if __name__=='__main__':render(Path(sys.argv[1]),Path(sys.argv[2]))
