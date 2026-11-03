"""Capture native HXT stdout and exact mesher versions for six research inputs."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'02_人体骨骼晶格肌肉方案/输出/shoulder_design'
HERE=Path(__file__).resolve().parent
MESHER=ROOT/'mujoco_arm/tools/build_shoulder_tpu_volume_mesh.py'
code="import sys,runpy;sys.path.insert(0,"+repr(str(ROOT/'.deps'))+");sys.argv[0]='build_shoulder_tpu_volume_mesh.py';runpy.run_path("+repr(str(MESHER))+",run_name='__main__')"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
runs=[]
for branch,subfolder,filename in [
    ('s1_tpu_surface_global_revision','', 'candidate_assembly_obj64.json'),
    ('s1_tpu_edge_collapse_revision','epsilon_0.001','candidate_assembly.json'),
    ('s1_tpu_local_fairing_revision','move_0.02','candidate_assembly.json'),
]:
    for module in ['infraspinatus','subscapularis']:
        assembly=D/branch/module/subfolder/filename
        dest=assembly.parent/'volume_mesh'/module/'preserve_h0.8_facetangle1e-06_algorithm10'
        dest.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(MESHER,dest/('generator_snapshot_'+sha(MESHER)[:12]+'.py'))
        # The previous uncaptured result is retained before replacing it.
        previous=dest/'mesh_attempt.json'
        if previous.exists():
            backup=dest/'before_native_capture'
            if backup.exists():raise RuntimeError('Existing capture run; inspect instead of restarting')
            backup.mkdir()
            for name in ['mesh_attempt.json','gmsh_log.txt']:
                p=dest/name
                if p.exists():shutil.copyfile(p,backup/name)
        command=[sys.executable,'-c',code,'--module',module,'--mode','preserve','--algorithm','10',
            '--facet-angle-tol','0.000001','--assembly',str(assembly)]
        print('Capturing',branch,module,flush=True)
        native_log=dest/'engine_stdout.txt'
        with native_log.open('wb') as stream:
            process=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        report=json.loads(previous.read_text(encoding='utf8'))
        row=dict(module=module,branch=branch,assembly_source=assembly.relative_to(ROOT).as_posix(),
            assembly_sha256=sha(assembly),attempt=previous.relative_to(ROOT).as_posix(),attempt_sha256=sha(previous),
            engine_log=native_log.relative_to(ROOT).as_posix(),engine_log_sha256=sha(native_log),
            process_exit_code=process.returncode,geometry_pass=report['geometry_tetrahedralization_pass'],
            numerical_quality_pass=report['numerical_elasticity_mesh_quality_pass'])
        runs.append(row)
        print(module,'exit',process.returncode,'geometry',row['geometry_pass'],'quality',row['numerical_quality_pass'],flush=True)
        (HERE/'captured_attempts.json').write_text(json.dumps(dict(runs=runs,complete=len(runs)==6,
            capture_script_sha256=sha(Path(__file__)),stage_pass=False),ensure_ascii=False,indent=2)+'\n',encoding='utf8')
