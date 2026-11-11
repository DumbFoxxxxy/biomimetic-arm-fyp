"""Reimported new material fit plus staged sternum-frame mounting."""
from itertools import combinations
import numpy as np
from build_s1_teres_origin_mount import A,ROOT,read,save,load_component
from build_channel_first_inventory import provenance


def main():
    folder=ROOT/'02_人体骨骼晶格肌肉方案/输出/channel_first/thorax_reaction_mount'
    file=folder/'assembly.json';assembly=read(file);parent=read(ROOT/assembly['parent']['file'])
    solids={k:load_component(v) for k,v in assembly['sources'].items()}
    key='thorax_reference:FJ3178';old_sternum=load_component(parent['sources'][key])
    new_keys=set(solids)-set(parent['sources']);changed={key,*new_keys};hits=[];inherited=[];pairs=0
    for a,b in combinations(solids,2):
        if a not in changed and b not in changed:continue
        pairs+=1;v=float(A.I.intersection_volume(solids[a],solids[b]))
        if v<=.001:continue
        other=b if a==key else a
        if key in (a,b) and other in parent['sources']:
            before=float(A.I.intersection_volume(old_sternum,solids[other]));extra=float(A.I.intersection_volume(solids[key]-old_sternum,solids[other]))
            if before>.001 and extra<=.001:
                inherited.append(dict(other=other,before_mm3=before,current_mm3=v,new_material_mm3=extra));continue
        hits.append(dict(left=a,right=b,overlap_mm3=v))
    save(folder/'export_screen.json',dict(assembly=provenance(file),tested_changed_pairs=pairs,
        new_interferences=hits,inherited_anatomical_junctions=inherited,
        generator=provenance(ROOT/'tools/audit_channel_first_thorax_mount.py'),
        scope='Exact reloaded new/modified rigid material. Original14 anatomical junction intersections are individually reported and new material separately tested. They remain unresolved whole-chest physical assembly; no whole-pair removal.'))
    assert not hits,hits
    fixed={k:v for k,v in solids.items() if k in parent['sources']};stages=[];normal=np.array([0.,-1.,0.])
    def stage(label,keys,axis,distance):
        failure=None;checked=0
        for delta in np.arange(distance,-.01,-.5):
            collisions=[]
            for k in keys:
                moving=solids[k].translate(axis*delta)
                for other,s in fixed.items():
                    v=float(A.I.intersection_volume(moving,s))
                    if v>.001:collisions.append(dict(moving=k,target=other,overlap_mm3=v))
            checked+=1
            if collisions:failure=dict(offset_mm=float(delta),hits=collisions);break
        stages.append(dict(stage=label,moving_keys=keys,fixed_keys=list(fixed),axis_source=axis.tolist(),
            initial_offset_mm=distance,sampling_mm=.5,checked_states=checked,first_interference=failure))
        if not failure:fixed.update({k:solids[k] for k in keys})
        print(label,'states',checked,'first',failure,flush=True);return failure is None
    success=True
    for i in (0,1):
        if not stage(f'spacer_{i}',[f'thorax_mount:spacer_{i}'],normal,30.):success=False;break
    module=['thorax_mount:reaction_frame']+[k for k in solids if k.startswith('thorax_port:')]
    if success:success=stage('preassembled_two_port_frame',module,normal,60.)
    if success:
        for i in (0,1):
            front=[f'thorax_mount:{suffix}_{i}' for suffix in ('bolt','front_washer')]
            rear=[f'thorax_mount:{suffix}_{i}' for suffix in ('nut','rear_washer')]
            if not stage(f'M3_{i}_front',front,normal,65.) or not stage(f'M3_{i}_rear',rear,-normal,45.):success=False;break
    save(folder/'installation_screen.json',dict(assembly=provenance(file),stages=stages,
        generator=provenance(ROOT/'tools/audit_channel_first_thorax_mount.py'),all_selected_stages_clear=success,
        scope='Neutral torso held in source pose, spacers manually held, port module already preassembled outside chest. Tests only mounting module/spacers/M3 groups. Does not validate port preassembly, crimping, threads/torque/tools, hand support,53-part anatomical chest joints or tabletop base.'))


if __name__=='__main__':main()
