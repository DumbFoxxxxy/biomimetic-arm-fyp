from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'.deps'),str(ROOT/'mujoco_arm/tools')]
from build_shoulder_tendon_segments import open_ended_tube,rod


class TendonPortTests(unittest.TestCase):
    def test_wire_passes_through_straight_end_caps(self):
        s=open_ended_tube([[0,0,0],[0,0,12]],1.8)
        self.assertEqual(len(s.decompose()),1)
        wire=rod([0,0,-3],[0,0,15],.3)
        self.assertLess(abs((s^wire).volume()),1e-6)
        self.assertGreater((s^rod([1.2,0,2],[1.2,0,10],.2)).volume(),.1)

    def test_bent_port_is_open_along_both_local_tangents(self):
        p=np.array([[120,-50,1300],[120,-50,1305],[124,-50,1309]],float)
        s=open_ended_tube(p,1.8)
        for a,b in [(p[0]-[0,0,3],p[1]),(p[1],p[-1]+np.array([1,0,1])*3)]:
            self.assertLess(abs((s^rod(a,b,.3)).volume()),1e-5)

    def test_no_wall_is_rejected(self):
        with self.assertRaises(ValueError):open_ended_tube([[0,0,0],[0,0,12]],.5,.85)


if __name__=='__main__':unittest.main()
