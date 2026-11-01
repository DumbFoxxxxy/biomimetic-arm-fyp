"""Attachment candidate regions of the anatomical arm: ONE loader.

`attachment_regions/<muscle>_<origin|insertion>_<bone>.npz` holds proximity
samples of the bone surface next to an anatomical muscle end, with per-sample
muscle face areas and distances.  Three studies already need the same
"representative attachment point" out of them - the force-line sign screen, the
humerus-side layout study and the muscle envelope study - so the reduction lives
here and nowhere else.

The reduction is deliberately simple and is NOT an anatomical claim:

    keep the samples at or below the median distance  (the near half)
    centroid = area-weighted mean of those bone-surface points

A candidate region is a proximity sample of a body-parts model, not a confirmed
attachment site (dissection or tendon footprint).  Every report that uses this
must say so.
"""
from pathlib import Path

import numpy as np

from frames import OUTPUT_DIR

REGION_DIR = OUTPUT_DIR / 'attachment_regions'
KINDS = ('origin', 'insertion')
BONES = ('scapula', 'clavicle', 'humerus')


def region_path(muscle_id, kind, bone):
    return REGION_DIR / ('%s_%s_%s.npz' % (muscle_id, kind, bone))


def load_region(muscle_id, kind, bone=None, near_half=True):
    """Load one candidate region.

    bone=None picks the first bone that has a file for this muscle/kind.
    Returns dict(points, weights, tris, bone, distance_min/max, area, n_samples)
    or None when no file exists.
    """
    if kind not in KINDS:
        raise ValueError('kind must be one of %r' % (KINDS,))
    bones = [bone] if bone else list(BONES)
    for name in bones:
        path = region_path(muscle_id, kind, name)
        if not path.is_file():
            continue
        with np.load(path) as d:
            points = d['bone_points_source_mm']
            weights = d['muscle_face_areas_mm2']
            distances = d['distances_mm']
            tris = d['bone_triangle_indices'] if 'bone_triangle_indices' in d else None
        keep = distances <= np.median(distances) if near_half else np.ones(len(points), bool)
        return {'points': points[keep], 'weights': weights[keep],
                'tris': None if tris is None else tris[keep],
                'bone': name, 'n_samples': int(keep.sum()),
                'n_all': int(len(points)),
                'distance_min_mm': float(distances[keep].min()),
                'distance_max_mm': float(distances[keep].max()),
                'area_mm2': float(weights[keep].sum())}
    return None


def centroid(region):
    """Area-weighted bone-surface centroid [mm] and its total area [mm^2]."""
    w = region['weights']
    return (region['points'] * w[:, None]).sum(0) / w.sum(), float(w.sum())
