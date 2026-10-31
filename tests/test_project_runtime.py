import os
import tempfile
import unittest
from pathlib import Path
from project_runtime import WorkspacePath,PROJECT_ROOT,project_path,load_mujoco_model,load_mujoco_xml
from project_runtime.paths import ROUTES

class PathTests(unittest.TestCase):
    def test_archived_routes(self):
        for old,new in ROUTES.items():
            self.assertEqual(str(project_path(old)),str(PROJECT_ROOT/new))
            self.assertEqual(str(project_path('D:/XJTLU/PSP/'+old)),str(PROJECT_ROOT/new))
    def test_join_in_pieces(self):
        self.assertEqual(str(WorkspacePath(PROJECT_ROOT)/'mujoco_arm'/'exports'/'anatomical_concept'),str(PROJECT_ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept'))
        self.assertEqual(str(WorkspacePath(PROJECT_ROOT).joinpath('mujoco_arm','exports','anatomical_concept')),str(PROJECT_ROOT/'02_人体骨骼晶格肌肉方案/输出/anatomical_concept'))
    def test_external_location(self):
        self.assertEqual(str(WorkspacePath('D:/XJTLU/isa_BP3D_4.0_obj_99')),str(Path('D:/XJTLU/isa_BP3D_4.0_obj_99')))
    def test_canonical_relative_location(self):
        self.assertEqual((project_path('骨骼')/'humerus.obj').relative_to(PROJECT_ROOT).as_posix(),'02_人体骨骼晶格肌肉方案/输入/骨骼/humerus.obj')
    def test_no_compatibility_junctions(self):
        import json
        for row in json.loads((PROJECT_ROOT/'tools/workspace_links.json').read_text(encoding='utf8')):
            self.assertFalse((PROJECT_ROOT/row['path']).exists(),row['path'])

class MuJoCoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (PROJECT_ROOT/'tmp').mkdir(exist_ok=True)
    def test_unicode_include_from_other_cwd(self):
        with tempfile.TemporaryDirectory(prefix='测试场景_',dir=PROJECT_ROOT/'tmp') as temp:
            p=Path(temp);(p/'资源').mkdir();(p/'子').mkdir()
            (p/'资源/四面体.obj').write_text('v 0 0 0\nv 0.01 0 0\nv 0 0.01 0\nv 0 0 0.01\nf 1 3 2\nf 1 2 4\nf 1 4 3\nf 2 3 4\n',encoding='utf8')
            (p/'子/包含.xml').write_text('<mujocoinclude><asset><mesh name="m" file="../资源/四面体.obj"/></asset><worldbody><geom type="mesh" mesh="m"/></worldbody></mujocoinclude>',encoding='utf8')
            scene=p/'入口.xml';scene.write_text('<mujoco><include file="子/包含.xml"/></mujoco>',encoding='utf8')
            cwd=os.getcwd()
            try:os.chdir(PROJECT_ROOT/'tools');model=load_mujoco_model(scene)
            finally:os.chdir(cwd)
            self.assertEqual(model.ngeom,1)
    def test_mesh_directory(self):
        with tempfile.TemporaryDirectory(prefix='场景_',dir=PROJECT_ROOT/'tmp') as temp:
            p=Path(temp);(p/'网格').mkdir()
            (p/'网格/a.obj').write_text('v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 3 2\nf 1 2 4\nf 1 4 3\nf 2 3 4\n',encoding='utf8')
            model=load_mujoco_xml('<mujoco><compiler meshdir="网格"/><asset><mesh name="m" file="a.obj"/></asset><worldbody><geom type="mesh" mesh="m"/></worldbody></mujoco>',base_dir=p)
            self.assertEqual(model.nmesh,1)
    def test_missing_asset(self):
        with self.assertRaises(FileNotFoundError):load_mujoco_xml('<mujoco><asset><mesh name="m" file="not_present.obj"/></asset></mujoco>')

if __name__=='__main__':unittest.main()
