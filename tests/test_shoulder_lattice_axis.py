"""Physical grid-direction checks, not a TPU material or actuation test."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'.deps'),str(ROOT/'mujoco_arm/tools')]
import build_shoulder_tpu as T
import mesh_geom


class LatticeAxisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.envelope=T.mf().Manifold.cube([200,200,200]).translate([-100,-100,-100])
        v,f=mesh_geom.manifold_mesh_to_arrays(cls.envelope.to_mesh64())
        cls.tri=v[f]
        cls.path=np.array([[0.,0.,0.],[0.,0.,60.]])

    def grid(self,**kwargs):
        return T.structured_lattice(self.envelope,self.tri,self.path,
            cell=(15.,10.,30.),i_range=(0,1),j_range=(0,1),**kwargs)

    def test_y_is_ten_mm_pitch_and_z_is_transverse(self):
        y,count_y,_=self.grid(path_axis='y')
        z,count_z,_=self.grid(path_axis='z')
        self.assertEqual(count_y,6)
        self.assertEqual(count_z,2)
        by=np.array(y.bounding_box());bz=np.array(z.bounding_box())
        self.assertGreater(by[3]-by[0],2.4*(bz[3]-bz[0]))
        self.assertGreater(y.volume(),0)
        self.assertGreater(z.volume(),0)

    def test_existing_call_keeps_z_mapping(self):
        default,nd,sd=self.grid()
        explicit,ne,se=self.grid(path_axis='z')
        self.assertEqual((nd,sd),(ne,se))
        self.assertLess(abs(default.volume()-explicit.volume()),1e-8)
        self.assertLess((default-explicit).volume()+(explicit-default).volume(),1e-8)


if __name__=='__main__':unittest.main()
