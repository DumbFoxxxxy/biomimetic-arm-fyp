"""M7 changed rigid installation paths and positive key seating checks.

Bare adapters fit before muscle modules. Keyed guide insertion is additionally
tested with all neutral muscle pieces present. This is not complete arm assembly.
"""
from build_guide_ground_frames import *

def run():
    a=json.loads((DEST7/'assembly.json').read_text());rows={r['key']:r for r in a['instances']};world={k:load(ROOT/r['file'])[0].translate(r['translation_display_mm']) for k,r in rows.items()};selected=json.loads((DEST/'paired_placement/selection.json').read_text())['placements'];paths=[];seats=[]
    input_paths=[DEST7/'assembly.json',DEST7/'design.json']+[ROOT/r['file'] for r in rows.values()];source_hashes={str(p):digest(p) for p in input_paths}
    def probe(label,moving,fixed,direction,length=25,base_offset=None):
        first=None;soft=[];samples=0;max_hard=0.
        for offset in np.r_[np.arange(length,5,-1),np.linspace(5,0,51)]:
            h=[];samples+=1
            for k in moving:
                s=world[k].translate(np.array(direction)*offset+(np.zeros(3) if base_offset is None else base_offset))
                for hit in hits(s,{j:world[j] for j in fixed}):
                    item={'moving':k,**hit,'offset_mm':float(offset)}
                    if rows[k]['kind']=='TPU' or rows[hit['part']]['kind']=='TPU':soft.append(item)
                    else:h.append(item);max_hard=max(max_hard,hit['volume_mm3'])
            if h and first is None:first={'offset_mm':float(offset),'hits':h}
        row={'label':label,'moving':moving,'direction':np.asarray(direction).tolist(),'samples':samples,'final_5mm_step_mm':.1,'other_step_mm':1,'first_rigid_obstruction':first,'rigid_clear':first is None,'soft_contact_flags':soft,'max_rigid_overlap_mm3':max_hard};paths.append(row);print(label,'clear',row['rigid_clear'],'soft',len(soft),flush=True);return row
    skeleton=[k for k in world if not k.startswith(('FJ1486_','FJ1477_'))];fitted=skeleton.copy()
    placements={p['id']:p for p in a['placements']}
    for fid in ['FJ1486','FJ1477']:
        n=np.array(placements[fid]['bone_anchors'][0]['bolt_axis']);adapter=fid+'_bone_adapter';normal=probe(fid+'_bare_extended_adapter_to_bone_before_muscle',[adapter],fitted,n,30)
        if SWEPT_ROOT and not normal['rigid_clear']:
            paths.pop();approaches=[normal];final=probe(fid+'_adapter_final_normal_5mm',[adapter],fitted,n,5);paths.pop()
            if final['rigid_clear']:
                ex=unit(np.array([1.,0,0])-n*n[0])
                for label,direction in [('medial',ex),('lateral',-ex)]:
                    lateral=probe(fid+'_adapter_side_staging_'+label,[adapter],fitted,direction,40,base_offset=n*5);paths.pop();approaches.append(lateral)
                    if lateral['rigid_clear'] and not lateral['soft_contact_flags']:
                        chosen={**lateral,'label':fid+'_bare_adapter_side_then_normal','legs':[lateral.copy(),final],'translation_staging_mm':(n*5).tolist(),'tested_approaches':[{k:r[k] for k in ['label','rigid_clear','first_rigid_obstruction']} for r in approaches]};paths.append(chosen);break
                else:paths.append(normal)
            else:normal['final_normal_leg']=final;paths.append(normal)
        fitted.append(adapter)
        for i in [0,1]:
            prefix=fid+'_bone_'+str(i);pair=[prefix+'_bolt',prefix+'_front_washer'];probe(prefix+'_front_bolt',pair,fitted,n,30);fitted+=pair
            pair=[prefix+'_nut',prefix+'_back_washer'];probe(prefix+'_back_nut',pair,fitted,-n,30);fitted+=pair
    for fid,s in selected.items():
        F=np.array(s['frame_columns']);Q=np.array(s['centre_display_mm']);keyparts=[fid+'_'+k for k in json.loads((TEMPLATE/'design.json').read_text())['parts']];keyhardware=[fid+'_'+k for k in local_receiver()[1]]
        # Full actual neutral surrounding context, no whole bone pair disabled.
        fixed=[k for k in world if k not in keyparts+keyhardware]
        probe(fid+'_preassembled_cassette_into_keyed_receiver',keyparts,fixed,F[:,1],25);fixed+=keyparts
        if LOW_ROOT:
            side=1 if fid=='FJ1486' else -1
            pair=[fid+'_key_nut',fid+'_key_front_washer'];options=[]
            for label,direction in [('outboard',F[:,0]*side),('minus_localY',-F[:,1]),('opposite_outboard',-F[:,0]*side),('plus_localY',F[:,1]),('diagonal_minusY',unit(F[:,0]*side-F[:,1])),('diagonal_plusY',unit(F[:,0]*side+F[:,1]))]:
                row=probe(fid+'_key_front_loose_nut_'+label,pair,fixed,direction,20);paths.pop();options.append(row)
                if row['rigid_clear'] and not row['soft_contact_flags']:break
            chosen=min(options,key=lambda r:(not r['rigid_clear'],len(r['soft_contact_flags']),r['max_rigid_overlap_mm3']));chosen['tested_directions']=[{k:r[k] for k in ['label','direction','rigid_clear','first_rigid_obstruction']} for r in options];paths.append(chosen);fixed+=pair
            pair=[fid+'_key_M3x18',fid+'_key_back_washer'];probe(fid+'_key_REAR_bolt',pair,fixed,-F[:,2],25)
        else:
            pair=[fid+'_key_nut',fid+'_key_back_washer'];probe(fid+'_key_rear_nut_first',pair,fixed,-F[:,2],20);fixed+=pair
            pair=[fid+'_key_M3x20',fid+'_key_front_washer'];probe(fid+'_key_front_bolt',pair,fixed,F[:,2],25)
        # Actual key rear seat is a positive stop, not a floating gap that a
        # tightened bolt would close. Patch lies in the frame, meets the fork.
        patch=(rod([0,-13,-5.2],[0,-13,-5],2.8,32)-rod([0,-13,-6],[0,-13,-4],1.7,32)).transform(np.c_[F,Q]);frame=world[fid+'_bone_adapter'];fork=world[fid+'_cassette_fork'];missing=max(0.,float((patch-frame).volume()));gap=float(patch.min_gap(fork,1));closing=intersection(patch,fork.translate(-F[:,2]*.01));seats.append({'muscle':fid,'missing_seat_patch_mm3':missing,'seat_patch_to_key_gap_mm':gap,'key_0_01mm_closing_overlap_mm3':closing,'pass':missing<.001 and gap<.0001 and closing>.01})
    result={'scope':'16 changed rigid installation paths with existing humeral bone holes; bare frame installed before muscles, guide insertion checked against neutral modules. No full muscle installation, tools, flexible handling, threads/preload/strength proof','paths':paths,'all_rigid_paths_clear':all(r['rigid_clear'] for r in paths),'TPU_contact_flag_count':sum(len(r['soft_contact_flags']) for r in paths),'positive_key_seat_checks':seats,'all_seats_pass':all(r['pass'] for r in seats),'source_sha256':source_hashes,'full_arm_assembly_verified':False,'frozen':False}
    assert all(digest(__import__('pathlib').Path(p))==h for p,h in source_hashes.items())
    (DEST7/'assembly_checks.json').write_text(json.dumps(result,indent=2),encoding='utf8');print('key seat checks',seats,flush=True)

if __name__=='__main__':run()
