"""Geometric bone-clear corridor screening, NOT a compliant muscle simulation.

Nodes are free numerical variables, not physical guide mounts or rail joints.
Only cap-pad endpoints are prescribed. A 2.1mm radius corridor is checked
against actual bone solids; other hardware and the full lattice remain outside
this scope. Never infer a working muscle or validated wire route from success.
"""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from build_direct_muscle_yoke import *
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from project_runtime import WorkspacePath as Path

class BoneDistance:
    def __init__(self,solid):
        self.solid=solid;mesh=solid.to_mesh64();v=np.asarray(mesh.vert_properties)[:,:3];self.tri=v[np.asarray(mesh.tri_verts)]
        self.tree=cKDTree(self.tri.mean(1));self.outside=np.array(solid.bounding_box())[3:]+[100,123,157]
    def query(self,points):
        _,ids=self.tree.query(points,k=min(128,len(self.tri)));tri=self.tri[ids];p=points[:,None,:];a,b,c=tri[:,:,0],tri[:,:,1],tri[:,:,2]
        ab=b-a;ac=c-a;n=np.cross(ab,ac);n2=(n*n).sum(2)
        projection=p-n*(((p-a)*n).sum(2)/np.maximum(n2,1e-30))[:,:,None]
        v0=ab;v1=ac;v2=projection-a;d00=(v0*v0).sum(2);d01=(v0*v1).sum(2);d11=(v1*v1).sum(2);d20=(v2*v0).sum(2);d21=(v2*v1).sum(2);den=np.maximum(d00*d11-d01*d01,1e-30)
        u=(d11*d20-d01*d21)/den;v=(d00*d21-d01*d20)/den;interior=(u>=0)&(v>=0)&(u+v<=1)
        candidates=[projection];valid=[interior]
        for x,y in [(a,b),(b,c),(c,a)]:
            edge=y-x;t=np.clip(((p-x)*edge).sum(2)/np.maximum((edge*edge).sum(2),1e-30),0,1)
            candidates.append(x+edge*t[:,:,None]);valid.append(np.ones_like(interior))
        nearest=np.stack(candidates,2);dist2=((p[:,:,None,:]-nearest)**2).sum(3);dist2=np.where(np.stack(valid,2),dist2,np.inf);flat=dist2.reshape(len(points),-1);index=flat.argmin(1);closest=nearest.reshape(len(points),-1,3)[np.arange(len(points)),index];length=np.sqrt(flat.min(1));sign=[]
        for point in points:
            # Segment ends beyond all bounds. Parity is authoritative solid
            # inside/outside; the nearby-triangle distance is an optimizer aid.
            sign.append(-1 if len(self.solid.ray_cast(point,self.outside))%2 else 1)
        sign=np.array(sign);gradient=(points-closest)/np.maximum(length[:,None],1e-12)*sign[:,None]
        return length*sign,gradient

def run():
    assert EXTENSOR_REANCHOR
    a=json.loads((DEST/'assembly.json').read_text());bykey={r['key']:r for r in a['instances']};definition=json.loads((BASE/'design.json').read_text());E=np.array(definition['joint_centre_display_mm']);axis=np.array(definition['joint_axis']);boundary=next(r for r in a['muscle_end_boundaries'] if r['muscle']=='FJ1477');A=np.array(boundary['proximal_pad_centre_display_mm']);B=np.array(boundary['distal_pad_zero_centre_display_mm']);offset=np.array(next(r['translation_display_mm'] for r in a['placements'] if r['id']=='FJ1477'));d=json.loads((muscle_folder('FJ1477')/'design.json').read_text());rest=np.array([A,*[np.array(g['centre_display_mm'])+offset for g in d['actual_short_tube_guides']],B]);N=len(rest);fractions=np.r_[0,np.cumsum(np.linalg.norm(np.diff(rest,axis=0),axis=1))];fractions/=fractions[-1]
    ports='--ports' in sys.argv;fixed_count=2 if ports else 1
    if ports:
        z=np.array(d['basis_columns_display'])[:,2]
        rest=np.array([A-z*5.15,A-z*2.15,*[np.array(g['centre_display_mm'])+offset for g in d['actual_short_tube_guides']],B+z*2.15,B+z*4.6]);N=len(rest);fractions=np.r_[0,np.cumsum(np.linalg.norm(np.diff(rest,axis=0),axis=1))];fractions/=fractions[-1]
    solids={key:load(ROOT/bykey[key]['file'])[0].translate(bykey[key]['translation_display_mm']) for key in ['FJ3368','FJ3391']};distance={key:BoneDistance(s) for key,s in solids.items()};inputs=[DEST/'assembly.json',muscle_folder('FJ1477')/'design.json']+[ROOT/bykey[k]['file'] for k in solids];hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};records=[];last=None;radius=2.1
    dense='--dense' in sys.argv
    angles=[0,90,120,131,132,135,140,145,150] if '--quick' in sys.argv else range(151)
    weights=list(np.eye(N))
    for i in range(N-1):
        for t in (np.arange(1,8)/8 if dense else [.5]):
            w=np.zeros(N);w[i]=1-t;w[i+1]=t;weights.append(w)
    weights=np.array(weights)
    for angle in angles:
        R=Rotation.from_rotvec(axis*np.radians(angle)).as_matrix();D=E+R@(B-E);end_targets=E+(rest[-fixed_count:]-E)@R.T;initial=rest+(D-B)*fractions[:,None]
        if last is not None:initial=last+(end_targets[-1]-last[-1])*fractions[:,None]
        initial[:fixed_count]=rest[:fixed_count];initial[-fixed_count:]=end_targets
        def unpack(x):return np.vstack([rest[:fixed_count],x.reshape(-1,3),end_targets])
        def objective(x):
            P=unpack(x);links=np.diff(P,axis=0);length=np.sqrt((links*links).sum(1)+1e-12);curvature=np.diff(P,n=2,axis=0);cost=float(length.sum()+.015*(curvature*curvature).sum());g=np.zeros_like(P);u=links/length[:,None];g[:-1]-=u;g[1:]+=u;g[:-2]+=.03*curvature;g[1:-1]-=.06*curvature;g[2:]+=.03*curvature
            return cost,g[fixed_count:-fixed_count].reshape(-1)
        # Midpoint samples supplement nodes during optimization. Final
        # acceptance is independent exact solid/capsule Boolean, not this SDF.
        def constraints(x):
            P=unpack(x);Q=weights@P;values=[];jacobians=[]
            for key,query in distance.items():
                local=(Q-E)@R+E if bykey[key]['owner']=='elbow' else Q
                value,g=query.query(local)
                if bykey[key]['owner']=='elbow':g=g@R.T
                J=weights[:,:,None]*g[:,None,:]
                values.extend(value-radius-(.3 if dense else .15));jacobians.append(J[:,fixed_count:-fixed_count].reshape(len(Q),-1))
            return np.array(values),np.vstack(jacobians)
        result=minimize(lambda x:objective(x)[0],initial[fixed_count:-fixed_count].reshape(-1),jac=lambda x:objective(x)[1],constraints=[{'type':'ineq','fun':lambda x:constraints(x)[0],'jac':lambda x:constraints(x)[1]}],method='SLSQP',options={'maxiter':180,'ftol':1e-7})
        P=unpack(result.x);pieces=[rod(p,q,radius,24) for p,q in zip(P,P[1:])]+[mf.Manifold.sphere(radius,16).translate(p) for p in P[1:-1]];corridor=mf.Manifold.batch_boolean(pieces,mf.OpType.Add);hits=[]
        for key,s in solids.items():
            target=s.transform(np.c_[R,E-R@E]) if bykey[key]['owner']=='elbow' else s;v=intersection(corridor,target)
            if v>.001:hits.append({'bone':key,'volume_mm3':v})
        record={'elbow_deg':int(angle),'optimizer_success':bool(result.success),'optimizer_message':str(result.message),'iterations':int(result.nit),'optimizer_constraint_min_mm':float(constraints(result.x)[0].min()),'centreline_display_mm':P.tolist(),'polyline_length_mm':float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum()),'radius_mm':radius,'exact_corridor_bone_intersections':hits,'bone_corridor_geometrically_clear':not hits,'actual_lattice_deformation_verified':False};records.append(record);last=P if not hits else None
        print('bone corridor',angle,'optimizer',result.success,'exactclear',not hits,'length',record['polyline_length_mm'],flush=True)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    stem='bone_corridor_quick' if '--quick' in sys.argv else 'bone_corridor_all_angles'
    report={'scope':'bone-only numerical centreline corridor feasibility; no new physical guide or external reaction point; NOT compliant equilibrium, actual guide tube pose, metal clearance, lattice deformation, actual wire route, friction or actuation verification','quick': '--quick' in sys.argv,'dense':dense,'ports':ports,'radius_mm':radius,'distance_method':'nearest128triangle centres plus exact ray parity for optimization; final independent manifold capsule-solid Booleans','endpoint_boundary':'two real rigid bore-axis points at each end prescribed' if ports else 'actual cap-pad centres prescribed, interior nodes entirely free; no pose constraints added to actual model','source_sha256':hashes,'poses':records,'frozen':False};(DEST/(stem+('_dense' if dense else '')+('_ports' if ports else '')+'.json')).write_text(json.dumps(report,indent=2),encoding='utf8')

if __name__=='__main__':run()
