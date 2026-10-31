"""Separate M3 rotor pocket experiment; never overwrite the M2 baseline."""

# Locate the project package from this script, independent of working directory.
import sys as _project_sys
from pathlib import Path as _NativeProjectPath
_project_root = next(p for p in _NativeProjectPath(__file__).resolve().parents if (p / 'project_runtime').is_dir())
if str(_project_root) not in _project_sys.path:
    _project_sys.path.insert(0, str(_project_root))
from project_runtime import load_mujoco_model as _load_project_model, load_mujoco_xml as _load_project_xml
from project_runtime import WorkspacePath as Path
import sys,json,hashlib,os,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'mujoco_arm/exports/anatomical_concept/elbow_muscles/mounted_elbow_candidate'
DEST=BASE/'distal_rotor_space/pocket_candidate'

def build():
    assembly=json.loads((BASE/'assembly.json').read_text(encoding='utf8'))
    substitutions={'FJ3368':'humerus_rotor_pocket','FJ3391':'forearm_with_rotor'}
    for r in assembly['instances']:
        if r['key'] in substitutions:
            r['file']=str((DEST/(substitutions[r['key']]+'.obj')).relative_to(ROOT));r['translation_display_mm']=[0,0,0]
    assembly['scope']='M3 explicit humeral rotor pocket and fused forearm rotor; provisional, no tendon yet'
    assembly['generated_components']={name:1 for name in substitutions.values()}
    assembly['source_M2_assembly_sha256']=hashlib.sha256((BASE/'assembly.json').read_bytes()).hexdigest()
    (DEST/'assembly.json').write_text(json.dumps(assembly,indent=2),encoding='utf8')
    for filename in ['arm.xml','neutral_view.xml']:
        tree=ET.parse(BASE/filename);root=tree.getroot();root.set('model','M3_rotor_pocket_experiment')
        for mesh in root.findall('./asset/mesh'):
            key=mesh.get('name');source=DEST/(substitutions[key]+'.obj') if key in substitutions else (BASE/mesh.get('file')).resolve()
            mesh.set('file',os.path.relpath(source,DEST).replace('\\','/'))
        ET.indent(root);tree.write(DEST/filename,encoding='utf-8',xml_declaration=True)
    scope=json.loads((BASE/'rigid_motion_scope.json').read_text(encoding='utf8'))
    scope['scope']='M3 candidate known rigid parts; actual fused rotor, no output tendon or compliant motion'
    (DEST/'rigid_motion_scope.json').write_text(json.dumps(scope,indent=2),encoding='utf8')
    print('Built separate pocket experiment, original M2 unchanged')

def run():
    if '--cad' in sys.argv:
        raise RuntimeError('Use export_rotor_pocket_native.py --humerus-only, then --canonical-mesh, import_rotor_cad_mesh.py and export_rotor_pocket_native.py; final Boolean microtriangle sewing is rejected.')
    if '--audit' in sys.argv:
        sys.argv+=['--muscle-mounts','--only-probes','elbow_flexion']
        import audit_joint_model as audit
        audit.DEST=DEST;audit.run();return
    if '--neutral' in sys.argv:
        import check_elbow_muscle_mounts as check
        check.DEST=DEST;check.run();return
    if '--exports' in sys.argv:
        import check_mount_exports as check
        check.DEST=DEST;check.run();return
    if '--mount-assembly' in sys.argv:
        import check_muscle_mount_assembly as check
        check.DEST=DEST;check.run();return
    build()

if __name__=='__main__':run()
