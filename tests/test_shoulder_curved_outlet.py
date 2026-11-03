"""Analytic curve/face checks independent of the anatomical search result."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mujoco_arm/tools'))
from design_shoulder_origin_outlet import curved_outlet


class CurvedOutletTests(unittest.TestCase):
    def test_circular_bend_and_exact_face_crossing(self):
        n=np.array([0.,0.,1.]);u=np.array([1.,0.,0.])
        result=curved_outlet(np.zeros(3),n,u,1.,60,3.)
        arc=result['arc'];centre=np.array([3.,0.,-np.sqrt(3)])
        np.testing.assert_allclose(np.linalg.norm(arc-centre,axis=1),3.,atol=1e-12)
        self.assertAlmostEqual(result['inner'][2],0.)
        self.assertAlmostEqual(result['outer'][2],1.)
        # Both ends have the intended analytic tangent: no ideal sharp corner.
        radius0=(arc[0]-centre)/3;radius1=(arc[-1]-centre)/3
        self.assertAlmostEqual(float(radius0@n),0.)
        self.assertAlmostEqual(float(radius1@result['direction']),0.)
        self.assertLessEqual(np.max(np.degrees(np.arccos(np.clip(
            np.sum((arc[:-1]-centre)*(arc[1:]-centre),axis=1)/9,-1,1)))),1.+1e-9)

    def test_rejects_bend_exceeding_available_inlet(self):
        self.assertIsNone(curved_outlet(np.zeros(3),np.array([0.,0.,1.]),
                                       np.array([1.,0.,0.]),1.,60,5.))


if __name__=='__main__':unittest.main()
