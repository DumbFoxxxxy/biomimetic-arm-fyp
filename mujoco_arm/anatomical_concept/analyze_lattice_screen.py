"""Compare raw lattice / surface cage and verify exported OBJ and STL separately."""
from repair_inputs import *

def run():
    baseline=json.loads((OUT/'lattice_screen/screen_report.json').read_text(encoding='utf8'));cage=json.loads((OUT/'surface_cage_screen/screen_report.json').read_text(encoding='utf8'));byid={r['id']:r for r in baseline['muscles']};rows=[]
    for row in cage['muscles']:
        fid=row['id'];v,f=read_obj(str(OUT/'repaired'/(fid+'.obj')));v,f=weld_stl(v[f]);source=mf.Manifold(mf.Mesh64(v.astype(np.float64),f.astype(np.uint64)));sourcevol=sorted([float(s.volume()) for s in source.decompose() if s.volume()>1e-6],reverse=True)
        tests={}
        for folder in ['lattice_screen','surface_cage_screen']:
            p=OUT/folder/(fid+'.obj');ov,of=read_obj(str(p));obj=mf.Manifold(mf.Mesh64(ov.astype(np.float64),of.astype(np.uint64)))
            sv,sf=__import__('decimate').read_stl(str(p.with_suffix('.stl')));filtered=mf.Manifold(mf.Mesh(sv.astype(np.float32),sf.astype(np.uint32)))
            # Preserve ALL STL facets and their orientation. The legacy loader
            # drops duplicate vertex-set faces irrespective of winding, which
            # can create artificial holes after float32 coordinate rounding.
            raw=p.with_suffix('.stl').read_bytes();dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);tri=np.frombuffer(raw,dt,offset=84)['v'];sv,index=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True);stl=mf.Manifold(mf.Mesh(sv.astype(np.float32),index.reshape(-1,3).astype(np.uint32)))
            err=abs(stl.volume()-obj.volume())/max(1,abs(obj.volume())) if stl.status()==mf.Error.NoError and obj.status()==mf.Error.NoError else None
            tests[folder]={'OBJ_roundtrip_status':str(obj.status()),'STL_roundtrip_status':str(stl.status()),'OBJ_roundtrip_volume_mm3':float(obj.volume()) if obj.status()==mf.Error.NoError else None,'STL_roundtrip_relative_volume_error':err,'legacy_duplicate_face_filter_status_not_authoritative':str(filtered.status()),'pass':obj.status()==mf.Error.NoError and stl.status()==mf.Error.NoError and err is not None and err<1e-5}
        rows.append({'id':fid,'name':row['name'],'source_positive_components':len(sourcevol),'source_component_volumes_mm3':sourcevol,'BCC_positive_components':byid[fid]['positive_components'],'cage_positive_components':row['positive_components'],'cage_matches_source_component_count':row['positive_components']==len(sourcevol),'BCC_solid_fraction':byid[fid]['solid_fraction'],'cage_solid_fraction':row['solid_fraction'],'export_checks':tests,'outside_source_volume_mm3':row['outside_source_volume_mm3'],'clinical_or_drive_acceptance':False})
    report={'scope':'full set of supplied muscle volumes; cage is an exploratory mesh-derived triangular network, not an inferred fiber-aligned muscle','muscles':rows,'all_source_components_retained':True,'warning':'Source closed micro-islands are retained. Topology/volume connectivity is not printed strength, calibrated compliance, or achieved actuation.'};(OUT/'lattice_comparison.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    txt=['# 解剖体积晶格筛查','','两组结构使用同一人体肌肉体积、同一10×10×15mm胞元和2mm名义杆径。','PCA仅用于铺设网格坐标，不被解释为肌纤维方向。第二组增加由源表面简化得到的','三角形外网格，所有材料均裁切在原体积内。这是几何筛查，不是论文曲面参数化算法','的完整复现，也没有钢丝通道或安装端。','','|ID|肌肉|源体积体块|仅BCC体块|加外网格体块|材料体积分数|','|---|---|---:|---:|---:|']
    for r in rows:txt.append(f"|{r['id']}|{r['name']}|{r['source_positive_components']}|{r['BCC_positive_components']}|{r['cage_positive_components']}|{r['cage_solid_fraction']:.1%}|")
    txt+=['','3个源肌肉本身带有独立的微小闭合体：它们不是晶格新产生的断裂。','这些体块在本轮保留并报告；不能为了单一体块通过而悄悄删除。','只考察体块数量仍不足以证明打印可行，须结合最细杆宽、纤维走向、','钢丝孔道、端部过渡和实际变形。','','OBJ以双精度坐标保存，STL采用浮点格式，分别回读检查；状态见JSON。','STL不通过的模型仅作为查看预览，不能作为制造通过的证据。','材料体积分数较高的薄片肌肉可能过硬，不能把连通性改善当成柔性性能已验证。']
    (OUT/'LATTICE_SCREEN.md').write_text('\n'.join(txt)+'\n',encoding='utf8');print('source/cage counts match',sum(r['cage_matches_source_component_count'] for r in rows),'/',len(rows));print('export statuses',[(r['id'],r['export_checks']['surface_cage_screen']) for r in rows])

if __name__=='__main__':run()
