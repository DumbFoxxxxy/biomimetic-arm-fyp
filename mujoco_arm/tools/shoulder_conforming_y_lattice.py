"""Shared-node cuboid lattice mapped into measured anatomical cross sections.

Paper Y follows the centreline; X/Z describe the section. Nominal dimensions
select cell counts, while local pitches deform with the input volume. Every
actual station/frame/node is recorded; nominal dimensions are not claimed as
constant or as an experimentally calibrated material model.
"""
import numpy as np
import build_shoulder_tpu as T


def build(envelope,path,cell,seed,section,record):
    P=np.asarray(path,float);seg=np.linalg.norm(np.diff(P,axis=0),axis=1);arc=np.r_[0,np.cumsum(seg)]
    count=max(3,int(np.ceil(arc[-1]/cell[1]))+1)
    xspan=np.ptp(section['estimated_middle_X_intersections_mm']);zspan=np.ptp(section['estimated_middle_Z_intersections_mm'])
    nx=max(1,min(4,int(np.floor(xspan/cell[0]+.5))));nz=max(1,min(6,int(np.floor(zspan/cell[2]+.5))))
    stations=[];frames=[];profiles=[]
    mf=T.mf()
    for k,target in enumerate(np.linspace(0.,arc[-1],count)):
        i=min(max(0,int(np.searchsorted(arc,target,side='left'))-1),len(P)-2)
        tangent=(P[i+1]-P[i])/seg[i];point=P[i]+tangent*(target-arc[i])
        if k==0:
            x=np.asarray(seed)-tangent*float(np.asarray(seed)@tangent);x/=np.linalg.norm(x)
        else:
            previous=frames[-1][:,1];cross=np.cross(previous,tangent);size=np.linalg.norm(cross)
            if size>1e-10:
                axis=cross/size;angle=np.arccos(np.clip(previous@tangent,-1.,1.));K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0.]])
                R=np.eye(3)+np.sin(angle)*K+(1.-np.cos(angle))*(K@K);x=R@frames[-1][:,0]
            else:x=frames[-1][:,0]
            x-=tangent*float(x@tangent);x/=np.linalg.norm(x)
        z=np.cross(x,tangent);F=np.column_stack([x,tangent,z]);frames.append(F);stations.append(point)
        # Ray-cast the actual closed envelope, not a projected box. The
        # first positive exit bounds the star-shaped section about the wire.
        def extent(direction):
            hits=envelope.ray_cast(point,point+direction*100.)
            distances=sorted(float((np.array(h.position)-point)@direction) for h in hits
                             if float((np.array(h.position)-point)@direction)>.01)
            return distances[0] if distances else None
        e=[extent(x),extent(-x),extent(z),extent(-z)]
        if any(v is None for v in e):
            raise ValueError('No closed cross-section exit at station '+str(k))
        profiles.append(dict(centre_source_mm=point.tolist(),paper_basis_columns_source=F.tolist(),
                             X_positive_negative_extent_mm=e[:2],Z_positive_negative_extent_mm=e[2:]))
    nodes={};node_rows=[]
    for k,(point,F,profile) in enumerate(zip(stations,frames,profiles)):
        for i in range(nx+1):
            for j in range(nz+1):
                a=-1.+2.*i/nx;b=-1.+2.*j/nz
                ex=profile['X_positive_negative_extent_mm'][0 if a>=0 else 1]
                ez=profile['Z_positive_negative_extent_mm'][0 if b>=0 else 1]
                offset=F[:,0]*(a*ex)+F[:,2]*(b*ez)
                if np.linalg.norm(offset)<1e-9:q=point.copy();radius=0.
                else:
                    direction=offset/np.linalg.norm(offset)
                    hits=envelope.ray_cast(point,point+direction*100.)
                    ds=sorted(float((np.array(h.position)-point)@direction) for h in hits
                              if float((np.array(h.position)-point)@direction)>.01)
                    if not ds:raise ValueError('No volume-mapped node boundary')
                    # Preserve grid interior fractions; boundary corners map
                    # to the section and are inset to keep shared nodes solid.
                    radius=.88*ds[0]*max(abs(a),abs(b));q=point+direction*radius
                nodes[i,j,k]=q;node_rows.append(dict(index=[i,j,k],source_mm=q.tolist(),radial_distance_mm=radius))
    pieces=[];centres=[];cells=[];live_nodes=set()
    for k in range(count-1):
        for i in range(nx):
            for j in range(nz):
                corners=[(i+di,j+dj,k+dk) for di in [0,1] for dj in [0,1] for dk in [0,1]]
                c=np.mean([nodes[v] for v in corners],axis=0);centres.append(c.tolist());live_nodes.update(corners)
                for vertex in corners:
                    piece=T.rod(c,nodes[vertex],1.)
                    if piece is not None:pieces.append(piece)
                pieces.append(mf.Manifold.sphere(1.,16).translate(c))
                cells.append(dict(index=[i,j,k],centre_source_mm=c.tolist(),corner_indices=[list(v) for v in corners]))
    # An open diamond pattern over boundary quads protects the inner lattice
    # without adding straight longitudinal skin rails. It is a geometry
    # candidate; its effect on stiffness and folding requires measurements.
    perimeter=[(i,0) for i in range(nx+1)]+[(nx,j) for j in range(1,nz+1)]+[(i,nz) for i in range(nx-1,-1,-1)]+[(0,j) for j in range(nz-1,0,-1)]
    surface_faces=[];surface_centres=[]
    for k in range(count-1):
        for j,(a,b) in enumerate(zip(perimeter,perimeter[1:]+perimeter[:1])):
            corners=[(*a,k),(*b,k),(*a,k+1),(*b,k+1)]
            centre=np.mean([nodes[v] for v in corners],axis=0);surface_centres.append(centre.tolist())
            for vertex in corners:
                piece=T.rod(centre,nodes[vertex],.75)
                if piece is not None:pieces.append(piece)
            pieces.append(mf.Manifold.sphere(.75,16).translate(centre))
            surface_faces.append(dict(axial_slab=k,perimeter_edge=j,centre_source_mm=centre.tolist(),corner_indices=[list(v) for v in corners]))
    record.update(construction='shared-node cuboid mapped by section rays',path_axis='y',requested_cell_xyz_mm=list(cell),
                  dimensions_are_nominal=True,cell_count_X_Z=[nx,nz],path_length_mm=float(arc[-1]),station_spacing_mm=float(arc[-1]/(count-1)),
                  surface_node_radial_fraction=.88,stations=profiles,nodes=node_rows,cells=cells,
                  cell_centres_source_mm=centres,live_corner_nodes_source_mm=[nodes[v].tolist() for v in sorted(live_nodes)],
                  surface_pattern='open diamond faces on the transported boundary grid',surface_strut_diameter_mm=1.5,
                  surface_faces=surface_faces,surface_junction_centres_source_mm=surface_centres,
                  all_axial_slabs_have_emitted_cells=True)
    return T.union_all(pieces),len(cells),0
