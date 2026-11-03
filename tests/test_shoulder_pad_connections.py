"""Regression for separated bosses that passed the old watertight-only gate."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'.deps'),str(ROOT/'mujoco_arm/tools')]
import build_shoulder_ends as BE


class ShoulderPadConnections(unittest.TestCase):
    def test_flat_pad_boss_is_joined_and_has_correct_height(self):
        solid=BE.disc([0,0,0],[0,0,1],4,boss_mm=3)
        self.assertEqual(len(solid.decompose()),1)
        self.assertAlmostEqual(solid.bounding_box()[5],7,places=6)
        self.assertEqual(str(solid.status()),'Error.NoError')

    def test_conforming_pad_boss_remains_joined_after_bone_cut(self):
        bone=BE.mf().Manifold.cube([40,40,2]).translate([-20,-20,-1])
        solid=BE.conforming_disc([0,0,0],[0,0,1],4,11,{'bone':bone},'own',{},boss_mm=3)
        self.assertEqual(len(solid.decompose()),1)
        self.assertAlmostEqual(solid.bounding_box()[5],7,places=6)
        self.assertLess(abs((solid ^ bone).volume()),1e-6)

    def test_rotated_pad_keeps_continuity(self):
        n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
        solid=BE.disc([120,-55,1300],n,4,boss_mm=3)
        self.assertEqual(len(solid.decompose()),1)


if __name__=='__main__':unittest.main()
