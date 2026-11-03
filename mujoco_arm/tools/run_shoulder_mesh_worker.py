"""Run native meshing in an isolated child with time/log resource budgets.

Operational limits are not geometric acceptance gates. A killed, crashed or
unfinished worker never counts as a mesh pass. Logs and source versions stay
available even when native code bypasses Python exception/finally handling.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
ENTRIES={'volume':'build_shoulder_tpu_volume_mesh.py','remesh':'remesh_shoulder_tpu_analysis.py'}


def run(entry,output,arguments,timeout_seconds=120.,maximum_log_bytes=8*1024*1024):
    output=output.resolve()
    if not output.is_relative_to(ROOT):raise ValueError('Worker output must stay inside project')
    if timeout_seconds<=0 or maximum_log_bytes<=0:raise ValueError('Resource budgets must be positive')
    output.mkdir(parents=True,exist_ok=True)
    if (output/'process_result.json').exists():raise ValueError('Existing result: use a new output directory')
    script=ROOT/'mujoco_arm/tools'/ENTRIES[entry]
    with script.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    shutil.copyfile(script,output/'worker_source.py')
    code='import sys,runpy;sys.path.insert(0,'+repr(str(ROOT/'.deps'))+');sys.argv[0]='+repr(script.name)+';runpy.run_path('+repr(str(script))+',run_name="__main__")'
    command=[sys.executable,'-c',code]+list(arguments)
    record=dict(entry=entry,worker_script=script.relative_to(ROOT).as_posix(),worker_script_sha256=digest,
        arguments=list(arguments),timeout_seconds=timeout_seconds,maximum_log_bytes=maximum_log_bytes,
        operational_limits_not_geometric_tolerances=True,worker_started=False,worker_terminal=False,
        adopted_as_current=False,stage_pass=False)
    log=output/'engine_stdout.txt';started=time.monotonic();reason=None;process=None
    try:
        with log.open('wb') as stream:
            process=subprocess.Popen(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
            record.update(worker_started=True,worker_pid=process.pid)
            while process.poll() is None:
                if log.stat().st_size>=maximum_log_bytes:reason='native_stdout_resource_limit'
                elif time.monotonic()-started>=timeout_seconds:reason='native_runtime_resource_limit'
                if reason:
                    process.kill();break
                time.sleep(.05)
            record['process_exit_code']=process.wait()
    except BaseException as error:
        reason='supervisor_interrupted_or_failed'
        record.update(supervisor_error_type=type(error).__name__,supervisor_error=str(error))
        if process is not None and process.poll() is None:process.kill()
        if process is not None:record['process_exit_code']=process.wait()
        raise
    finally:
        if process is not None:record['worker_terminal']=process.poll() is not None
        record.update(stop_reason=reason,elapsed_s=time.monotonic()-started,
            worker_report_may_be_incomplete=reason is not None or record.get('process_exit_code')!=0)
        if log.exists():
            with log.open('rb') as stream:record['engine_stdout_sha256']=hashlib.file_digest(stream,'sha256').hexdigest()
            record['engine_stdout_bytes']=log.stat().st_size
        (output/'process_result.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('Native worker exit',record.get('process_exit_code'),'stop reason',reason,flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--entry',choices=ENTRIES,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout-seconds',type=float,default=120.)
    parser.add_argument('--maximum-log-bytes',type=int,default=8*1024*1024)
    parser.add_argument('worker_arguments',nargs=argparse.REMAINDER)
    args=parser.parse_args();forwarded=args.worker_arguments
    if forwarded[:1]==['--']:forwarded=forwarded[1:]
    result=run(args.entry,args.output,forwarded,args.timeout_seconds,args.maximum_log_bytes)
    if result['stop_reason'] or result['process_exit_code']!=0:raise SystemExit(3)
