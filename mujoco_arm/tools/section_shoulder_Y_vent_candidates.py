"""Export exact planar solid sections for source/vent comparison diagrams."""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_shoulder_retained_joint_fit as A

DEST=A.DESIGN/'s1_y_muscle_pair/vent_experiment'

def run():
    file=DEST/'vent_trials.json';d=json.loads(file.read_text(encoding='utf8'))
    independent_file=DEST/'independent_checks.json';checked=json.loads(independent_file.read_text(encoding='utf8'))
    if checked['trial_report_sha256']!=A.hashed(file):raise ValueError('Stale verification')
    sections=[]
    for row in d['results']:
        if row['module']!='infraspinatus' or not row['all_cavities_opened_in_search']:continue
        source=row['parent_TPU'];old,_,_=A.load_checked(A.ROOT/source['file'],source['translation_source_mm'],source['sha256'])
        mesh=row['mesh'];actual,_,_=A.load_checked(A.ROOT/mesh['file'],mesh['translation_source_mm'],mesh['sha256'])
        chosen=row['steps'][0]['selected'];p,q=np.array(chosen['path_source_mm'])
        x=A.D.unit(q-p)
        void=next(s for s in old.decompose() if s.volume() < -1.)
        points=np.array(void.to_mesh64().vert_properties)[:,:3]
        y=points.mean(axis=0)-p;y-=x*(x@y);y=A.D.unit(y);z=np.cross(x,y)
        F=np.column_stack([x,y,z]);T=np.c_[F.T,-F.T@p]
        local=(points-p)@F
        low=local[:,:2].min(axis=0)-2.;high=local[:,:2].max(axis=0)+2.
        # Include both endpoints so the cut corridor is visible.
        high[0]=max(high[0],float(np.linalg.norm(q-p))+1.)
        box=A.I._manifold().Manifold.cube([high[0]-low[0],high[1]-low[1],40.]).translate([low[0],low[1],-20.])
        panels=[]
        for title,solid in [('Source: closed pocket',old),('Trial: local opening',actual)]:
            section=(solid.transform(T)^box).slice(0.)
            panels.append(dict(label=title,contours_local_mm=[np.asarray(poly).tolist() for poly in section.to_polygons()]))
        sections.append(dict(aperture_diameter_candidate_mm=row['aperture_diameter_candidate_mm'],
            source_sha256=source['sha256'],trial_sha256=mesh['sha256'],
            source_to_section_affine=T.tolist(),cut_segment_local_mm=[[0.,0.],[float(np.linalg.norm(q-p)),0.]],
            bounds_local_mm=[low.tolist(),high.tolist()],panels=panels,
            interpretation='exact zero-height planar slice; independent trial rows use their own selected cut planes'))
    report=dict(trial_report_sha256=A.hashed(file),independent_checks_sha256=A.hashed(independent_file),
        script_sha256=A.hashed(Path(__file__).resolve()),sections=sections,
        not_verified='planar diagram is not a printing, strength, channel deformation or driving result')
    (DEST/'section_geometry.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
