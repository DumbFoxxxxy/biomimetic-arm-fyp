"""Provide actual insertion clearance for the fused rotor, keep fork ears."""
from check_rotor_pocket import ROOT,BASE,DEST
from build_elbow_muscle_mounts import load,rod,mf,np
from build_SC_connection import write
import json,hashlib,shutil
def run():
    selection=json.loads((DEST/'selection.json').read_text(encoding='utf8'))
    if (DEST/'rotor_assembly_relief.json').exists():raise RuntimeError('Assembly relief already recorded; rebuild source candidate before reapplying')
    elbow=ROOT/'mujoco_arm/exports/anatomical_concept/joint_model/SC_connection/AC_compact/GH_connection/elbow_connection'
    prior=json.loads((elbow/'assembly_and_retention_checks.json').read_text(encoding='utf8'))
    frame=json.loads((DEST.parent/'space_screen.json').read_text(encoding='utf8'));E=np.array(frame['centre_display_mm']);axis=np.array(frame['axis']);direction=np.array(prior['chosen_cartridge_path']['direction'])
    geometry=selection['geometry'];station=geometry['axis_station_mm'];r=geometry['nominal_steel_centre_radius_mm'];lo,hi=station-2.35,station+2.35
    start,end=sorted([np.sign(station)*7.6,station]);cuts=[]
    for distance in range(46):
        origin=E+direction*distance
        cuts+=[rod(origin+axis*lo,origin+axis*hi,r+1.35),rod(origin+axis*(start-.35),origin+axis*(end+.35),8.35)]
    hum=load(DEST/'humerus_rotor_pocket.obj')[0];revised=(hum-mf.Manifold.batch_boolean(cuts,mf.OpType.Add)).as_original().simplify(.01)
    components=[c for c in revised.decompose() if c.volume()>.001]
    if revised.status()!=mf.Error.NoError or len(components)!=1:
        rejected={'candidate':'swept straight rotor insertion recess','status':str(revised.status()),'positive_component_volumes_mm3':[float(c.volume()) for c in components],'component_bounds_mm':[list(c.bounding_box()) for c in components],'accepted':False,'source_geometry_changed':False,'reason':'the required straight insertion opening breaks humeral connectivity; no fragment deleted, no assembly/motion pass claimed','next_design':'replace integral rotor by actual detachable torque-coupled rotor, or independently examine other install directions'}
        (DEST/'rejected_insertion_relief.json').write_text(json.dumps(rejected,indent=2),encoding='utf8');print(rejected)
        return
    assert max(abs(lo),abs(hi),abs(start-.35),abs(end+.35))<34
    archive=DEST/'before_rotor_insertion_relief';archive.mkdir(exist_ok=True)
    for name in ['humerus_rotor_pocket.obj','humerus_rotor_pocket.stl','cad_input.json','joint_assembly_checks.json','bone_motion_audit.json','rotor_cad_roundtrip.json']:
        if (DEST/name).exists():shutil.copy2(DEST/name,archive/name)
    part=write(DEST/'humerus_rotor_pocket.stl',revised);assert part['components']==1
    data=json.loads((DEST/'cad_input.json').read_text(encoding='utf8'));data['parts']=[part if p['name']==part['name'] else p for p in data['parts']]
    (DEST/'cad_input.json').write_text(json.dumps(data),encoding='utf8')
    report={'reason':'bare forearm with fused rotor collided during approach at28mm despite clear final pose; actual localized insertion recess added','source_humerus_OBJ_sha256':hashlib.sha256((archive/'humerus_rotor_pocket.obj').read_bytes()).hexdigest(),'translation_direction':direction.tolist(),'construction_translation_range_mm':[0,45],'construction_step_mm':1,'radial_and_axial_nominal_clearance_mm':.35,'axis_extent_mm':[min(lo,start-.35),max(hi,end+.35)],'fork_ears_begin_at_abs_axis_mm':34,'ear_axial_ranges_untouched':True,'positive_material_components':1,'additional_humeral_removal_mm3':float((hum-revised).volume()),'full_motion_and_assembly_recheck_required':True,'frozen':False}
    (DEST/'rotor_assembly_relief.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(report)
if __name__=='__main__':run()
