from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'mujoco_arm/tools'))
from analyze_shoulder_local_actuation import tension_moments,positive_spanning
import interference as I


class LocalActuationTests(unittest.TestCase):
    def test_positive_and_negative_cartesian_lines_span(self):
        m=np.c_[np.eye(3),-np.eye(3)]
        r=positive_spanning(m)
        self.assertTrue(r['strict_positive_span'])
        self.assertAlmostEqual(r['normalized_positive_margin'],1/6,places=6)

    def test_positive_only_directions_do_not_span(self):
        result=positive_spanning(np.eye(3))
        self.assertFalse(result['strict_positive_span'])
        self.assertTrue(all(x>0 for x in result['strict_separating_plane']['all_column_projections_mm']))

    def test_balanced_but_planar_directions_do_not_control_three_axes(self):
        m=np.array([[1,-1,0,0],[0,0,1,-1],[0,0,0,0]])
        self.assertFalse(positive_spanning(m)['strict_positive_span'])

    def test_virtual_work_sign_and_radian_units(self):
        A=np.array([[1.,0.,0.]]);B=np.array([[0.,1.,0.]]);G=np.zeros(3);k=np.array([0.,0.,1.])
        M,L=tension_moments(A,B,G)
        self.assertAlmostEqual(float(k@M[:,0]),-1/np.sqrt(2))
        for deg in [0,30,80]:
            T=I.rotation_about(G,k,deg);moving=B@T[:,:3].T+T[:,3]
            M,_=tension_moments(A,moving,G);h=1e-4
            def length(theta):
                T=I.rotation_about(G,k,theta)
                return np.linalg.norm(B@T[:,:3].T+T[:,3]-A)
            derivative=(length(deg+h)-length(deg-h))/(2*np.radians(h))
            self.assertAlmostEqual(float(k@M[:,0]),-derivative,places=7)


if __name__=='__main__':unittest.main()
