"""Losslessly archive a numerically rejected response without parsing huge pose arrays."""
import argparse
import gzip
import hashlib
import json
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    parser.add_argument('--result-name',default='short_response');args=parser.parse_args()
    assert args.result_name and all(c.isalnum() or c in '_-' for c in args.result_name)
    allowed=(ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_elbow_supply_dynamic').resolve()
    folder=(ROOT/args.folder).resolve();assert folder.is_relative_to(allowed)
    source=folder/(args.result_name+'.json');original=provenance(source)
    with source.open('rb') as stream:
        prefix=stream.read(65536);stream.seek(max(0,source.stat().st_size-65536));tail=stream.read()
    start=prefix.index(b'"trace": [');end=tail.rindex(b'\n  ],\n  "adaptive_rotation_trigger_deg"')
    metadata=json.loads(prefix[:start]+b'"trace": []'+tail[end+len(b'\n  ]'):])
    count=0;last=None;backward=0
    with source.open('rb') as stream:
        for line in stream:
            line=line.lstrip()
            if line.startswith(b'"time_s": '):
                time=float(line.split(b':',1)[1].strip().rstrip(b','))
                backward+=int(last is not None and time<=last);last=time;count+=1
    reasons=[]
    if metadata['warning_counts']:reasons.append('solver_warnings')
    if metadata['numerical_axial_guard_hit']:reasons.append('numerical_axial_guard')
    if metadata['maximum_recorded_rotation_increment_deg']>1.:reasons.append('insufficient_angular_sampling')
    if backward:reasons.append('nonmonotonic_time_after_solver_reset')
    if last is None or last+2*metadata['timestep_s']<metadata['duration_s']:reasons.append('requested_duration_not_completed')
    assert reasons,'Refuse to archive a viable response with this failure-only tool'
    target=source.with_suffix('.json.gz');assert not target.exists(),target
    with source.open('rb') as inp,gzip.open(target,'wb',compresslevel=6) as out:
        while chunk:=inp.read(1024*1024):out.write(chunk)
    digest=hashlib.sha256()
    with gzip.open(target,'rb') as inp:
        while chunk:=inp.read(1024*1024):digest.update(chunk)
    assert digest.hexdigest()==original['sha256'],'Lossless archive verification failed'
    snapshot=None;generator=ROOT/metadata['generator']['file']
    if provenance(generator)==metadata['generator']:
        diagnostics=allowed/'diagnostics';diagnostics.mkdir(exist_ok=True)
        file=diagnostics/(generator.stem+'_'+original['sha256'][:12]+'.py');file.write_bytes(generator.read_bytes())
        snapshot=provenance(file)
    record=dict(original_uncompressed=original,archive=provenance(target),original_bytes=source.stat().st_size,
        archive_bytes=target.stat().st_size,decompressed_sha256=digest.hexdigest(),rejection_reasons=reasons,
        metadata_without_poses=metadata,saved_states=count,last_saved_time_s=last,nonmonotonic_time_intervals=backward,
        original_generator_snapshot=snapshot,generator=provenance(ROOT/'tools/archive_channel_first_failed_response.py'),
        scope='Numerically invalid response, not a mechanical motion result. Every original byte retained in gzip and SHA256 verified before raw copy removal. Original code snapshot retained if still current.')
    save(folder/(args.result_name+'_numerical_rejection.json'),record)
    source.unlink()  # Sole redundant raw copy; exact verified gzip retained above.
    print('Rejected numerical response archived losslessly:',record['original_bytes'],'->',record['archive_bytes'],'bytes;',reasons,flush=True)


if __name__=='__main__':main()
