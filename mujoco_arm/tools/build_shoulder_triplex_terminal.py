"""Shared deltoid insertion prototype from the unified three-wire routes.

Preserves the R11 seat and M3 pitch. Exports a cap, TPU pad and three uncrimped
stop blanks, not a complete muscle or a proven moving/retaining assembly.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import design_shoulder_port_paths as D
import build_shoulder_ends as BE
import check_shoulder_assembly as C
import shoulder_pads
import mesh_geom
import interference as I
import build_shoulder_tendon_segments as T


def export_part(path,solid,origin):
    local=I.clean_solid(solid.translate(-origin))
    v,f=mesh_geom.manifold_mesh_to_arrays(local.to_mesh())
    mesh_geom.write_stl(str(path),v,f,header='S1 triplex terminal prototype')
    restored=C.load_stl(path);v,f=mesh_geom.read_stl(str(path));v,f=mesh_geom.weld(v,f)
    _,boundary,nm=mesh_geom.edge_census(f)
    error=abs(restored.volume()-local.volume())/local.volume()
    if len(restored.decompose())!=1 or boundary or nm or error>1e-5:
        raise ValueError('Invalid terminal export: '+str(path))
    return restored.translate(origin),dict(file=path.name,translation_source_mm=origin.tolist(),
        material_components=1,boundary_edges=boundary,nonmanifold_edges=nm,
        relative_volume_error=error,volume_mm3=float(restored.volume()),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    source=D.DEST/'routes_separated.json';design=json.loads(source.read_text(encoding='utf8'))
    routes=[r for r in design['routes'] if r['module']=='deltoid']
    seats=json.loads((D.DESIGN/'s1_end_pads.json').read_text(encoding='utf8'))['seats']
    seat=next(s for s in seats if len(s['heads'])==3)
    centre=np.array(seat['centre_source_mm']);n=D.unit(seat['normal_source']);basis=shoulder_pads.pad_basis(n)
    inner=np.array(routes[0]['ports']['distal']['inner_face_source_mm'])-np.array(routes[0]['ports']['distal']['offset_from_seat_mm'])
    dest=D.DEST/'triplex_terminal';dest.mkdir(exist_ok=True)
    original=D.frames.OUTPUT_DIR/'shoulder_parts/FJ1468_distal+FJ1513_distal+FJ1467_distal_pad.stl'
    parent=C.load_stl(original);pad=parent
    cap=BE.disc(inner,n,seat['cap_thickness_mm'],radius=seat['cap_radius_mm'],holes=False)
    for sign in [-1,1]:
        point=inner+basis[:,0]*sign*7
        cap=cap-T.rod(point-n*3,point+n*(seat['cap_thickness_mm']+3),1.7)
    lumens=[];wires={};stops={}
    for row in routes:
        head=row['head'];offset=np.array(row['ports']['distal']['offset_from_seat_mm'])
        pad_height=seat['lift_mm']+seat['pad_thickness_mm']
        # Local dense TPU collar, joined to the existing pad without a coplanar gap.
        collar=I._manifold().Manifold.cylinder(3.01,2.5,2.5,48)
        collar=collar.transform(np.c_[basis,centre+offset+n*(pad_height-.01)])
        pad=pad+collar
        full=np.array(row['full_wire_path_source_mm'])
        # Interface clearance drill: 1.72mm diameter (core lumen stays1.70).
        # A small explicit overcut removes knife-edge remnants where the new
        # curved passages meet the existing central opening in this prototype.
        lumen=D.outer_envelope(full,.86);lumens.append(lumen)
        wires[head]=D.outer_envelope(full,.3)
        out=np.array(row['ports']['distal']['outer_face_source_mm'])
        # Each wire has its own blind counterbore and nominal uncrimped sleeve.
        cap=cap-T.rod(out-n*2.6,out+n*.5,1.4)
        sleeve=T.rod(out-n*2.55,out-n*.05,1.25)-T.rod(out-n*3,out+n*.2,.45)
        stops[head]=sleeve
    lumen=I._manifold().Manifold.batch_boolean(lumens,I._manifold().OpType.Add)
    pad=pad-lumen;cap=cap-lumen
    result={};parts={}
    for name,s in [('TPU_pad',pad),('rigid_cap',cap)]+[(h+'_stop_blank',s) for h,s in stops.items()]:
        restored,row=export_part(dest/(name+'.stl'),s,centre);parts[name]=restored;result[name]=row
    checks={}
    for h,wire in wires.items():
        checks[h]={name:I.intersection_volume(wire,s) for name,s in parts.items()}
    bone_dir=D.frames.PROJECT_ROOT/'02_人体骨骼晶格肌肉方案/输入/骨骼'
    bone_hits={name:{key:I.intersection_volume(s,D.SP.load_solid(bone_dir/file)) for key,file in D.SP.BONES.items()}
               for name,s in parts.items()}
    rigid_soft_overlap=I.intersection_volume(parts['TPU_pad'],parts['rigid_cap'])
    stop_capture={}
    for head in stops:
        sleeve=parts[head+'_stop_blank']
        stop_capture[head]={str(distance):I.intersection_volume(
            sleeve.translate(-n*distance),parts['rigid_cap']) for distance in [0.,.04,.06]}
    clear=all(v<=.001 for row in checks.values() for v in row.values())
    print('Triplex5 parts exported; wire intersections',checks,flush=True)
    print('Bone intersections',bone_hits,'pad/cap overlap',rigid_soft_overlap,flush=True)
    report=dict(scope='static three-port shared insertion candidate, NOT full S1 or physical crimp validation',
        routes_source=str(source.relative_to(D.frames.PROJECT_ROOT)),
        routes_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        original_pad_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),parts=result,
        nominal_crimp_counterbore_radius_mm=1.4,nominal_crimp_depth_mm=2.6,
        interface_clearance_bore_radius_mm=.86,
        stop_blank_outer_diameter_mm=2.5,stop_blank_length_mm=2.5,wire_diameter_mm=.6,
        wire_intersections_mm3=checks,bone_intersections_mm3=bone_hits,pad_cap_overlap_mm3=rigid_soft_overlap,
        stop_inward_translation_intersections_mm3=stop_capture,
        all_nominal_wire_probes_clear=clear,
        parent_pad_material_removed_mm3=float((parent-pad).volume()),
        pad_material_added_mm3=float((pad-parent).volume()),
        stage_pass=False,not_verified=['source muscle/connection geometry not regenerated yet','neighbour modules and fixture',
            'crimp grip, fasteners, retention and strength','installation paths','deformation, friction and drive range'])
    (dest/'check.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print('Inward stop capture probes at0/.04/.06mm:',stop_capture,flush=True)


if __name__=='__main__':main()
