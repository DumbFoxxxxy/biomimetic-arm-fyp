from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mujoco_arm/tools'))
import compact_shoulder_terminal as CT
import build_shoulder_tendon_segments as T
import interference as I


class CompactTerminalTests(unittest.TestCase):
    def setUp(self):
        self.n=np.array([0.,0.,1.]);self.cap,self.stop,self.spec=CT.build([0,0,0],self.n,[1,0,0])

    def test_separate_wire_and_M3_paths_clear_actual_parts(self):
        wire=T.rod([0,0,-2],[0,0,5],.3)
        bolt=T.rod([5,0,-2],[5,0,5],1.5)
        self.assertLess(I.intersection_volume(wire,self.cap),1e-8)
        self.assertLess(I.intersection_volume(wire,self.stop),1e-8)
        self.assertLess(I.intersection_volume(bolt,self.cap),1e-8)
        self.assertLess(I.intersection_volume(bolt,self.stop),1e-8)

    def test_stop_has_clearance_then_positive_shoulder_capture(self):
        self.assertLess(I.intersection_volume(self.stop,self.cap),1e-8)
        self.assertLess(I.intersection_volume(self.stop.translate(-self.n*.04),self.cap),1e-8)
        self.assertGreater(I.intersection_volume(self.stop.translate(-self.n*.06),self.cap),.01)
        self.assertEqual(sum(p.volume()>1e-9 for p in self.cap.decompose()),1)


if __name__=='__main__':unittest.main()
