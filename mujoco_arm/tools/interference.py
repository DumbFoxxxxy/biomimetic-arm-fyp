"""Rigid interference and probe sweeps with ONE set of thresholds.

Before this module every study carried its own::

    def probe(label, moving, fixed, direction, length): ...   # 13 different copies
    if volume > .001: ...                                     # literal threshold

with the shared primitive living in one unrelated builder.  Those conventions are
now here and nowhere else.

What this module does NOT do (state it in any report you generate):

* no compliant deformation, friction, preload or assembly sequence;
* a sweep is discrete sampling plus a bisection on the boundary, not a proof;
* ``gap`` is capped at ``GAP_CAP_MM`` - "capped" means "at least that clear".

Requires manifold3d (see .cursor/skills/muscle-to-tpu/reference.md for the local
install note).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

#: penetration volume above which a pair is reported as a hit [mm3]
HARD_MM3 = 0.001
#: clearance queries are capped here; a bigger value means "at least this clear" [mm]
GAP_CAP_MM = 3.0
#: angular refinement of the first-contact boundary [deg]
BOUNDARY_DEG = 0.1

KINDS = ('rigid', 'tpu', 'wire')


def _manifold():
    try:
        import manifold3d
        return manifold3d
    except ImportError as exc:            # pragma: no cover - environment specific
        raise ImportError(
            'manifold3d is required for interference checks. '
            'Put the local .deps directory on PYTHONPATH; see '
            '.cursor/skills/muscle-to-tpu/reference.md section 3.1.'
        ) from exc


# --------------------------------------------------------------------- solids
@dataclass
class Solid:
    """A named piece of geometry taking part in a check.

    kind: 'rigid' -> any penetration is a hard hit
          'tpu'   -> recorded separately as a soft contact flag
          'wire'  -> same as 'tpu' (cable inside a compliant body)
    """

    key: str
    geom: object
    kind: str = 'rigid'

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError('kind must be one of %r, got %r' % (KINDS, self.kind))

    @property
    def is_soft(self):
        return self.kind in ('tpu', 'wire')


# ------------------------------------------------------------------ primitives
def intersection_volume(a, b):
    """Boolean intersection volume [mm3] with a bounding-box pre-reject."""
    aa = np.asarray(a.bounding_box(), dtype=float)
    bb = np.asarray(b.bounding_box(), dtype=float)
    if np.any(np.maximum(aa[:3] - bb[3:], bb[:3] - aa[3:]) >= 0.0):
        return 0.0
    return max(0.0, float((a ^ b).volume()))


def gap(a, b, cap=GAP_CAP_MM):
    """Minimum surface gap [mm], capped at `cap`."""
    return float(a.min_gap(b, cap))


def rotation_about(centre, axis, degrees):
    """3x4 affine rotating about (centre, axis) by `degrees`."""
    centre = np.asarray(centre, dtype=float)
    axis = np.asarray(axis, dtype=float)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        raise ValueError('zero-length axis')
    axis = axis / n
    t = math.radians(degrees)
    k = np.array([[0.0, -axis[2], axis[1]],
                  [axis[2], 0.0, -axis[0]],
                  [-axis[1], axis[0], 0.0]])
    R = np.eye(3) + math.sin(t) * k + (1.0 - math.cos(t)) * (k @ k)
    return np.c_[R, centre - R @ centre]


# ----------------------------------------------------------------- contacts
def contact(moving, fixed, hard_mm3=HARD_MM3):
    """Penetration of one solid against a list of solids.

    Returns dict(r_hits=[...], soft=[...], pass=bool) where `pass` means
    "no rigid penetration above the threshold".
    """
    rigid, soft = [], []
    for item in fixed:
        if item.key == getattr(moving, 'key', None):
            continue
        volume = intersection_volume(moving.geom, item.geom)
        if volume <= hard_mm3:
            continue
        record = {'part': item.key, 'kind': item.kind, 'volume_mm3': volume}
        (soft if item.is_soft else rigid).append(record)
    return {'rigid_hits': rigid, 'soft': soft, 'pass': not rigid}


def sweep(axis, centre, moving_geom, fixed: Sequence[Solid], angles_deg,
          *, moving_key='moving', hard_mm3=HARD_MM3, refine=True,
          boundary_deg=BOUNDARY_DEG, gap_cap=None):
    """Sweep `moving_geom` about (centre, axis) and record contacts per angle.

    `moving_geom` must be a manifold3d Manifold in the zero pose.
    `gap_cap` adds a per-angle minimum gap to each row when given.
    """
    angles = [float(a) for a in angles_deg]
    rows = []
    first_hit = None
    previous_clear = None
    worst = None

    def evaluate(angle):
        T = rotation_about(centre, axis, angle)
        moved = moving_geom.transform(T)
        res = contact(Solid(moving_key, moved), fixed, hard_mm3=hard_mm3)
        row = {'angle_deg': angle,
               'rigid_hits': res['rigid_hits'],
               'soft': res['soft'],
               'pass': res['pass']}
        if gap_cap is not None:
            row['min_gap_mm'] = min((gap(moved, s.geom, gap_cap) for s in fixed),
                                    default=float(gap_cap))
        return row

    for angle in angles:
        row = evaluate(angle)
        rows.append(row)
        total = sum(h['volume_mm3'] for h in row['rigid_hits'])
        if total > 0 and (worst is None or total > worst['volume_mm3']):
            worst = {'angle_deg': angle, 'volume_mm3': total}
        if not row['pass'] and first_hit is None:
            first_hit = angle
            if refine and previous_clear is not None:
                lo, hi = previous_clear, angle
                while abs(hi - lo) > boundary_deg:
                    mid = 0.5 * (lo + hi)
                    if evaluate(mid)['pass']:
                        lo = mid
                    else:
                        hi = mid
                first_hit = hi
        if row['pass']:
            previous_clear = angle

    return {
        'axis': list(map(float, np.asarray(axis, dtype=float))),
        'centre': list(map(float, np.asarray(centre, dtype=float))),
        'thresholds': {'hard_mm3': hard_mm3,
                       'boundary_deg': boundary_deg if refine else None,
                       'gap_cap_mm': gap_cap},
        'angles_deg': angles,
        'rows': rows,
        'first_hit_deg': first_hit,
        'worst': worst,
        'pass': first_hit is None,
        'scope': 'rigid geometry only; no compliance, friction, preload or assembly sequence',
        'boundary_is_sampled': True,
    }


# --------------------------------------------------------------------- export
def _clean(obj):
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    return obj


def write_report(path, report):
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(_clean(report), fh, ensure_ascii=False, indent=2)
    return path


# ------------------------------------------------------------------ self test
def _cube(size, offset=(0, 0, 0)):
    mf = _manifold()
    return mf.Manifold.cube(list(size)).translate(list(offset))


def self_test():
    ok = True

    def check(name, value, expected, tol):
        nonlocal ok
        good = abs(value - expected) <= tol
        print('[%s] %-34s %.6f (expect %.6f)' % ('PASS' if good else 'FAIL', name, value, expected))
        ok &= good

    # 1) two 10 mm cubes overlapping 5 mm in x -> 500 mm3
    a = _cube((10, 10, 10))
    b = _cube((10, 10, 10), (5, 0, 0))
    check('intersection 5 mm overlap', intersection_volume(a, b), 500.0, 1e-6)

    # 2) touching only -> 0
    c = _cube((10, 10, 10), (10, 0, 0))
    check('intersection touching', intersection_volume(a, c), 0.0, 1e-9)

    # 3) separated by 2 mm -> gap 2
    d = _cube((10, 10, 10), (12, 0, 0))
    check('gap 2 mm', gap(a, d), 2.0, 1e-6)

    # 4) gap cap
    e = _cube((10, 10, 10), (30, 0, 0))
    check('gap capped', gap(a, e, cap=3.0), 3.0, 1e-9)

    # 5) rotation helper: +90 deg about z maps (1,0,0) -> (0,1,0)
    T = rotation_about([0, 0, 0], [0, 0, 1], 90.0)
    p = T @ np.array([1.0, 0.0, 0.0, 1.0])
    check('rotation x->y (x)', p[0], 0.0, 1e-12)
    check('rotation x->y (y)', p[1], 1.0, 1e-12)

    # 6) sweep: thin rod (half width 2 mm, length 40 mm) about z; 1 mm obstacle
    #    centred at radius 20 mm, 30 deg.  Contact starts at
    #    30 - asin((2 + 0.5) / 20) ~= 22.82 deg.
    rod = _cube((40, 4, 4)).translate((0, -2, -2))       # x in [0,40], y,z in [-2,2]
    ang = math.radians(30.0)
    obs_c = np.array([20 * math.cos(ang), 20 * math.sin(ang), 0.0])
    obs = Solid('obstacle', _cube((1, 1, 1)).translate(obs_c - 0.5))
    expected = 30.0 - math.degrees(math.asin(2.5 / 20.0))

    def clear_at(angle):
        moved = rod.transform(rotation_about([0, 0, 0], [0, 0, 1], angle))
        return contact(Solid('rod', moved), [obs])['pass']

    rep = sweep([0, 0, 1], [0, 0, 0], rod, [obs], np.arange(0, 91, 1.0), refine=True)
    got = rep['first_hit_deg']
    good = got is not None and abs(got - expected) <= 1.0
    print('[%s] sweep first hit %.3f deg (expect %.3f +-1.0)'
          % ('PASS' if good else 'FAIL', got if got is not None else float('nan'), expected))
    ok &= good

    # 7) bisection: just before the reported boundary must be clear, just after must hit
    if got is not None:
        before = clear_at(got - 2 * BOUNDARY_DEG)
        after = clear_at(got + 2 * BOUNDARY_DEG)
        good = before and not after
        print('[%s] boundary bracketed: clear@-%.1f=%s  hit@+%.1f=%s'
              % ('PASS' if good else 'FAIL', 2 * BOUNDARY_DEG, before, 2 * BOUNDARY_DEG, after))
        ok &= good

    # 8) TPU parts are reported as soft, not as a pass-blocking hit
    tpu = Solid('muscle', _cube((10, 10, 10), (5, 0, 0)), kind='tpu')
    res = contact(Solid('rigid1', _cube((10, 10, 10))), [tpu])
    good = res['pass'] and len(res['soft']) == 1 and not res['rigid_hits']
    print('[%s] TPU overlap -> soft only: %s' % ('PASS' if good else 'FAIL', res))
    ok &= good

    print('\nRESULT:', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    import sys
    if '--self-test' in sys.argv or len(sys.argv) == 1:
        sys.exit(self_test())
    print(__doc__)
