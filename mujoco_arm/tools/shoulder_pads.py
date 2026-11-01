"""The bone-side seat envelope used by the S1 interface studies.

One definition, because three questions are answered with it: do seats overlap each
other, do they intrude into bone or bearing hardware, and how tall is the stack the
end cap has to clear.  Changing the pad here changes every answer.

The envelope is deliberately crude:

    a cylinder, radius 11.0 mm, thickness 4.0 mm, grown OUTWARD along the local
    bone-surface normal from the seat centre

It has no chamfer, no bracket, no bolt head, no tendon stub, and it does not know
about the 8.15 mm pad+cap stack that §5.7 puts on top of it.  Every report must say
so: a pass on this envelope is a pass on a flat disc, nothing more.

Coordinates are source BodyParts3D mm.  (The G1 joint parts are stored in display
mm; convert them with ``frames.SHIFT`` before mixing.)
"""
import math

import numpy as np

import frames

#: seat envelope [mm]
PAD_RADIUS_MM = 11.0
PAD_THICKNESS_MM = 4.0
#: seat fastener pattern [mm] - 2 x M3 at 14 mm pitch, from the elbow design
BOLT_PITCH_MM = 14.0
BOLT_HOLE_RADIUS_MM = 1.7
#: the holes must sit inside the seat
MIN_RADIUS_TWO_BOLT_MM = BOLT_PITCH_MM / 2.0 + BOLT_HOLE_RADIUS_MM


def pad_basis(normal):
    """Orthonormal basis with the third column along `normal`."""
    n = frames.unit(normal)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = frames.unit(np.cross(n, helper))
    y = np.cross(n, x)
    return np.stack([x, y, n], axis=1)


def build_pad(centre_source, normal_source, radius=PAD_RADIUS_MM,
              thickness=PAD_THICKNESS_MM, segments=48):
    """Cylinder seat sitting on the bone surface, growing outward along the normal.

    `Manifold.cylinder(h)` spans z from 0 to h, so the plate's inner face is already at the
    centre; translating by `thickness / 2` lifted every pad 2 mm off its seat (for a 4 mm plate).
    That offset is what made `lift_to_clear` report 0 on most seats - the pad was already 2 mm
    clear of the bone before any lift was applied.
    """
    import interference as I
    mf = I._manifold()
    R = pad_basis(normal_source)
    solid = mf.Manifold.cylinder(thickness, radius, radius, segments)
    return solid.transform(np.c_[R, np.asarray(centre_source, dtype=float)])


def pad_far_radius(centre_source, normal_source, centre_ref,
                   radius=PAD_RADIUS_MM, thickness=PAD_THICKNESS_MM):
    """Largest distance from `centre_ref` to any point of the seat envelope [mm].

    The seat centre can sit inside a sphere while its rim corners stick out, so
    anything that asks "does this clear a radius" must use this, not the centre.
    """
    c = np.asarray(centre_source, dtype=float)
    R = pad_basis(normal_source)
    worst = 0.0
    for t in (0.0, thickness):
        for k in range(72):
            a = 2.0 * math.pi * k / 72.0
            p = (c + R[:, 2] * t + R[:, 0] * radius * math.cos(a)
                 + R[:, 1] * radius * math.sin(a))
            worst = max(worst, float(np.linalg.norm(p - centre_ref)))
    return worst


def far_radius_about(centre_source, normal_source, centre_ref):
    return pad_far_radius(centre_source, normal_source, centre_ref)


def lift_to_clear(centre_source, normal_source, bone_solid, limit=8.0,
                  radius=PAD_RADIUS_MM, thickness=PAD_THICKNESS_MM):
    """Smallest outward offset along the normal with no bone intrusion [mm].

    A flat disc on a curved bone either floats on the high points or digs in at the
    rim; the lift says how much floating the flat form costs.  Returns
    (lift, residual_intrusion_at_lift); lift is None when even `limit` is not enough.
    """
    import interference as I

    def intrusion(offset):
        c = np.asarray(centre_source, float) + np.asarray(normal_source, float) * offset
        return I.intersection_volume(build_pad(c, normal_source, radius=radius,
                                               thickness=thickness), bone_solid)

    if intrusion(0.0) <= I.HARD_MM3:
        return 0.0, 0.0
    lo, hi = 0.0, limit
    if intrusion(hi) > I.HARD_MM3:
        return None, float(intrusion(hi))
    while hi - lo > 0.02:
        mid = 0.5 * (lo + hi)
        if intrusion(mid) <= I.HARD_MM3:
            hi = mid
        else:
            lo = mid
    return hi, 0.0
