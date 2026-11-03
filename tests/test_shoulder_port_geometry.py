from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mujoco_arm/tools'))
from design_shoulder_port_paths import trim_path,segment_distance,path_distance


class PortGeometryTests(unittest.TestCase):
    def test_trim_arclength_across_corner(self):
        p=trim_path([[0,0,0],[3,0,0],[3,4,0]],1,6)
        np.testing.assert_allclose(p,[[1,0,0],[3,0,0],[3,3,0]])

    def test_crossed_and_skew_segments(self):
        self.assertAlmostEqual(segment_distance([-1,0,0],[1,0,0],[0,-1,0],[0,1,0]),0)
        self.assertAlmostEqual(segment_distance([-1,0,0],[1,0,0],[0,-1,2],[0,1,2]),2)

    def test_parallel_and_degenerate_segments(self):
        self.assertAlmostEqual(segment_distance([0,0,0],[2,0,0],[1,1,0],[3,1,0]),1)
        self.assertAlmostEqual(segment_distance([1,0,0],[1,0,0],[0,0,0],[3,0,0]),0)

    def test_reversed_path_distance(self):
        p=np.array([[0,0,0],[3,0,0],[3,4,0]])
        self.assertAlmostEqual(path_distance(p,p[::-1]),0)

    def test_sum_endpoint_roundoff_but_not_real_overrun(self):
        p=np.random.default_rng(24).normal(size=(50,3)).cumsum(axis=0)
        total=np.linalg.norm(np.diff(p,axis=0),axis=1).sum()
        trimmed=trim_path(p,1.,total)
        np.testing.assert_allclose(trimmed[-1],p[-1],atol=1e-12)
        with self.assertRaises(ValueError):trim_path(p,1.,total+1e-5)


if __name__=='__main__':unittest.main()
