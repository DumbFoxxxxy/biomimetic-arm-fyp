"""Real-data regression for mujoco_arm/tools/interference.py.

Reproduces the three pin-capture witnesses recorded in KEEPER_REVIEW.md (M10-R1):

    原位及外移 0.09 mm 时不碰压片
    外移 0.11 mm 时被压片挡住，约 0.025 mm3
    反向移动 0.01 mm 时被原垫圈挡住，约 0.087 mm3

Run:
    set PYTHONPATH=<repo>\\.deps
    python mujoco_arm/tools/regression_interference.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import interference as I
import mesh_geom

BASE = (Path(__file__).resolve().parents[2] / '02_人体骨骼晶格肌肉方案' / '输出' /
        'anatomical_concept' / 'direct_muscle_interfaces' / 'mechanical_yoke_reanchored_pair')
KEEPERS = BASE / 'paired_swivel_support' / 'captive_keepers'
CAP = BASE / 'cap_swivel_flexor'

#: documented values, mm3
EXPECT_CLEAR_009 = 0.0
EXPECT_HIT_011 = 0.025
EXPECT_PULLOUT_001 = 0.087


def load(path):
    mf = I._manifold()
    V, F = mesh_geom.load_obj(str(path))
    return mf.Manifold(mf.Mesh(V.astype(np.float32), F.astype(np.uint32)))


def main():
    pin_path = KEEPERS / 'FJ1486_pin_1.obj'
    pin = load(pin_path)
    keeper = load(KEEPERS / 'FJ1486_keeper.obj')
    washer = load(CAP / 'washer_1.obj')

    V, _ = mesh_geom.load_obj(str(pin_path))
    centre = V.mean(axis=0)
    _, _, vt = np.linalg.svd(V - centre, full_matrices=False)
    axis = vt[0] / np.linalg.norm(vt[0])

    # The pin is backed by the cap washer on one side: that defines "inward".
    inward = None
    for sign in (1.0, -1.0):
        moved = pin.translate(list(axis * sign * 0.01))
        if I.intersection_volume(moved, washer) > 0.001:
            inward = axis * sign
            break
    if inward is None:
        print('FAIL: cannot determine the pin axis sign from the cap washer')
        return 1

    outward = -inward
    checks = []

    def witness(name, offset, target, expected, tol, must_be_zero=False):
        moved = pin.translate(list(outward * offset))
        volume = I.intersection_volume(moved, target)
        if must_be_zero:
            good = volume <= 0.001
        else:
            good = abs(volume - expected) <= tol
        checks.append({'name': name, 'offset_mm': offset, 'volume_mm3': volume,
                       'expected_mm3': expected if not must_be_zero else 0.0,
                       'pass': bool(good)})
        print('[%s] %-46s %8.4f mm3 (expect %s)'
              % ('PASS' if good else 'FAIL', name, volume,
                 '<= 0.001' if must_be_zero else '%.3f +- %.3f' % (expected, tol)))

    witness('outward +0.09 mm vs keeper (capture clear)', 0.09, keeper,
            EXPECT_CLEAR_009, 0.0, must_be_zero=True)
    witness('outward +0.11 mm vs keeper (blocked)', 0.11, keeper,
            EXPECT_HIT_011, 0.002)
    witness('inward 0.01 mm vs cap washer (pullout)', -0.01, washer,
            EXPECT_PULLOUT_001, 0.002)

    ok = all(c['pass'] for c in checks)
    report = {'scope': 'rigid capture witnesses only; no force, vibration or looseness claim',
              'pin_axis': axis.tolist(), 'inward_axis': inward.tolist(),
              'checks': checks, 'pass': ok}
    if '--json' in sys.argv:
        out = Path(sys.argv[sys.argv.index('--json') + 1])
        I.write_report(out, report)
        print('report written:', out)

    print('\nRESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
