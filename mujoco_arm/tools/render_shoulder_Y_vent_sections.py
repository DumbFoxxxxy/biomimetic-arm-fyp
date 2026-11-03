"""Render recorded solid contours with parity holes, without a CAD surrogate."""
import json
from pathlib import Path
from PIL import Image, ImageChops, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design/s1_y_muscle_pair/vent_experiment'

def run():
    d=json.loads((DEST/'section_geometry.json').read_text(encoding='utf8'))
    width=1200;row_height=290;top=120;height=top+len(d['sections'])*row_height+70
    img=Image.new('RGB',(width,height),'#f5f3ed');draw=ImageDraw.Draw(img)
    fontfile=Path('C:/Windows/Fonts/arial.ttf')
    font=ImageFont.truetype(str(fontfile),18) if fontfile.exists() else ImageFont.load_default()
    title=ImageFont.truetype(str(fontfile),25) if fontfile.exists() else font
    draw.text((28,20),'Infraspinatus: source and vent trial sections',font=title,fill='#163b43')
    draw.text((28,58),'Actual STL slices. Teal = TPU; pale = void/outside; orange = cut centreline.',font=font,fill='#364c52')
    draw.text((28,86),'Each row uses its own recorded cutting plane; compare left and right within the row.',font=font,fill='#364c52')
    for i,row in enumerate(d['sections']):
        low,high=row['bounds_local_mm'];dx=high[0]-low[0];dy=high[1]-low[1]
        scale=min(500/max(dx,1e-9),210/max(dy,1e-9))
        for j,panel in enumerate(row['panels']):
            left=30+j*590;yy=top+i*row_height
            draw.text((left,yy),f'{row["aperture_diameter_candidate_mm"]:g} mm | {panel["label"]}',font=font,fill='#163b43')
            ox=left+(550-dx*scale)/2.;oy=yy+35+(220-dy*scale)/2.
            def point(v):return (ox+(v[0]-low[0])*scale,oy+(high[1]-v[1])*scale)
            # Even/odd polygon parity preserves all holes regardless of
            # winding; no filling of enclosed shells occurs in the diagram.
            mask=Image.new('1',img.size,0)
            for contour in panel['contours_local_mm']:
                layer=Image.new('1',img.size,0);ImageDraw.Draw(layer).polygon([point(v) for v in contour],fill=1)
                mask=ImageChops.logical_xor(mask,layer)
            img.paste('#287b7e',(0,0,width,height),mask)
            draw=ImageDraw.Draw(img)
            for contour in panel['contours_local_mm']:
                pts=[point(v) for v in contour];draw.line(pts+[pts[0]],fill='#134e58',width=1)
            if j==1:
                draw.line([point(v) for v in row['cut_segment_local_mm']],fill='#d77c20',width=3)
            bx=left+12;by=yy+254;draw.line([(bx,by),(bx+5*scale,by)],fill='#163b43',width=3)
            draw.text((bx+5*scale+8,by-12),'5 mm',font=font,fill='#163b43')
    draw.text((28,height-48),'Geometry comparisons only. Default S1-Y1 unchanged; no manufacturing or drive pass.',font=font,fill='#364c52')
    img.save(DEST/'source_vs_vent_sections.png')

if __name__=='__main__':run()
