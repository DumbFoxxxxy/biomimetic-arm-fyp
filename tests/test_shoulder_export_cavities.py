from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mujoco_arm/tools'))
import interference as I
from build_shoulder_unified_muscle import discard_microscopic_positive_islands


class ExportCavityTests(unittest.TestCase):
    def test_removes_only_debris_and_keeps_enclosed_cavity(self):
        mf=I._manifold().Manifold
        body=mf.cube([10,10,10])-mf.cube([4,4,4]).translate([3,3,3])
        dust=mf.cube([.01,.01,.01]).translate([20,0,0])
        restored,removed=discard_microscopic_positive_islands(body+dust)
        self.assertEqual(len(removed),1)
        self.assertAlmostEqual(restored.volume(),936.,places=8)
        self.assertEqual(sum(p.volume()<0 for p in restored.decompose()),1)
        self.assertEqual(sum(p.volume()>0 for p in restored.decompose()),1)

    def test_real_disconnected_piece_still_fails_connectivity(self):
        mf=I._manifold().Manifold
        body=mf.cube([10,10,10])+mf.cube([1,1,1]).translate([20,0,0])
        restored,removed=discard_microscopic_positive_islands(body)
        self.assertEqual(removed,[])
        self.assertEqual(sum(p.volume()>0 for p in restored.decompose()),2)


if __name__=='__main__':unittest.main()
