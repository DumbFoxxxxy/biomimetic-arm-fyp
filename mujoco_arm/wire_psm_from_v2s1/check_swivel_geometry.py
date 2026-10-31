"""Independent geometric invariants for the finite-radius sheave calculation."""
import numpy as np
from common import swivel_route
def run():
    rng=np.random.default_rng(4);errors=[];length_errors=[];count=0
    for k in range(200):
        A=np.array([-100.,0,0]);G=np.zeros(3);B=rng.uniform(-150,150,3);r=swivel_route(A,G,B)
        if r is None:continue
        p=r['path'];radius=p[-2]-r['center'];out=p[-1]-p[-2];errors.extend([abs(np.linalg.norm(radius)-8),abs(out@radius)/max(np.linalg.norm(out),1)])
        length_errors.append(abs(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()-100-r['D_mm']));count+=1
    assert max(errors)<1e-8
    assert max(length_errors)<.01
    print({'samples':count,'max_radius_or_tangent_error_mm':max(errors),'max_polyline_arc_length_error_mm':max(length_errors)},flush=True)
if __name__=='__main__':run()
