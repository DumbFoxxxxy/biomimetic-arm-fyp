"""Projected actual open failure patches with detected edge/facet crossings."""
from pathlib import Path
import json
import numpy as np
from PIL import Image,ImageDraw,ImageFont
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y2_transition_stock/volume_mesh'

def run():
    data=json.loads((BASE/'diagnostics.json').read_text(encoding='utf8'))
    image=Image.new('RGB',(1200,750),'#f5f4ef');draw=ImageDraw.Draw(image)
    path=Path('C:/Windows/Fonts/arial.ttf');font=ImageFont.truetype(str(path),18) if path.exists() else ImageFont.load_default()
    title=ImageFont.truetype(str(path),25) if path.exists() else font
    draw.text((25,18),'Actual TPU boundary: localized volume-meshing failures',font=title,fill='#203c48')
    draw.text((25,58),'Open diagnostic patches; source anatomy/TPU unchanged. Markers indicate local edge/facet intersections.',font=font,fill='#203c48')
    # One fixed orthographic view for actual geometry, not a deformed surrogate.
    a=np.deg2rad(28.);b=np.deg2rad(35.)
    Rz=np.array([[np.cos(a),-np.sin(a),0.],[np.sin(a),np.cos(a),0.],[0.,0.,1.]])
    Rx=np.array([[1.,0.,0.],[0.,np.cos(b),-np.sin(b)],[0.,np.sin(b),np.cos(b)]])
    rotation=Rx@Rz
    for i,row in enumerate(data['localized_failures']):
        folder=(ROOT/row['open_patch']).parent;mesh=np.load(folder/'input_boundary.npz');V=mesh['nodes_local_mm'];F=mesh['triangles'];centre=np.array(row['engine_point_local_mm'])
        tri=V[F];low=tri.min(axis=1);high=tri.max(axis=1);r=row['localization_radius_mm']
        mask=np.all(low<=centre+r,axis=1)&np.all(high>=centre-r,axis=1)
        world=(tri[mask]-centre)@rotation.T
        xy=world[:,:,:2].reshape(-1,2);minimum=xy.min(axis=0);maximum=xy.max(axis=0)
        scale=min(510/max(maximum[0]-minimum[0],1e-12),420/max(maximum[1]-minimum[1],1e-12))
        left=25+i*590;top=150;ox=left+270-(minimum[0]+maximum[0])/2*scale;oy=top+220+(minimum[1]+maximum[1])/2*scale
        def project(p):
            q=(np.asarray(p)-centre)@rotation.T
            return (ox+q[0]*scale,oy-q[1]*scale)
        draw.text((left,110),row['module']+' | '+str(len(row['strict_nonadjacent_segment_facet_crossings']))+' local crossings',font=font,fill='#203c48')
        for j in np.argsort(world[:,:,2].mean(axis=1)):
            points=[(ox+p[0]*scale,oy-p[1]*scale) for p in world[j]]
            draw.polygon(points,fill='#b4ced2',outline='#527b88')
        for crossing in row['strict_nonadjacent_segment_facet_crossings']:
            edge=[project(p) for p in crossing['edge_points_local_mm']];draw.line(edge,fill='#c36822',width=3)
            x,y=project(crossing['intersection']['point_local_mm']);draw.ellipse((x-4,y-4,x+4,y+4),fill='#b33e30')
        y=640;draw.line([(left+15,y),(left+15+.1*scale,y)],fill='#203c48',width=3)
        draw.text((left+22+.1*scale,y-12),'0.1 mm',font=font,fill='#203c48')
        draw.text((left,680),'Projected source surface; not a stress or displacement plot.',font=font,fill='#203c48')
    image.save(BASE/'localized_boundary_failures.png')

if __name__=='__main__':run()
