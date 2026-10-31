"""131 actual neutral manufactured meshes, no invented muscle motion."""
from build_distal_swivels import *
def run():
    assembly=json.loads((BASE/'assembly.json').read_text(encoding='utf8'))
    replaced={fid+'_distal_rigid_cap' for fid in ['FJ1486','FJ1477']};assembly['instances']=[r for r in assembly['instances'] if r['key'] not in replaced]
    for fid in ['FJ1486','FJ1477']:
        folder=DEST/fid;design=json.loads((folder/'design.json').read_text(encoding='utf8'))
        for name in design['parts']:
            assembly['instances'].append({'key':fid+'_swivel_'+name,'file':str((folder/(name+'.obj')).relative_to(ROOT)),'translation_display_mm':[0,0,0],'owner':'distal_unvalidated','kind':'rigid','source_module':fid,'m5_module':fid,'m5_part':name})
    assert len(assembly['instances'])==131
    assembly.update(scope='M5 neutral packaging with real spherical end terminals; no whole-muscle deformation/drive/installed path validation',replaced_parent_keys=sorted(replaced),parent_assembly_sha256=hashlib.sha256((BASE/'assembly.json').read_bytes()).hexdigest(),frozen=False)
    (DEST/'assembly.json').write_text(json.dumps(assembly,indent=2),encoding='utf8');print('M5 neutral scene instances',len(assembly['instances']))
if __name__=='__main__':run()
