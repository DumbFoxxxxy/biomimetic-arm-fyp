"""Joint frames and coordinate transforms for the anatomical arm.

Single source of truth, so studies stop re-deriving (and mis-reading) these.

Coordinate chain
    source BodyParts3D  --(-SHIFT)-->  display      SHIFT = [-140, -70, 1100]
    display             --(-E) @ F --> muscle-local (u, v, h)

Muscle-local basis (columns of F)
    z = unit(E_elbow - G_gh)        humerus axis, GH -> elbow
    x = unit(cross(z, [0, 1, 0]))
    y = cross(z, x)

Verified against FJ1486_profile.json: max|F - F_profile| = 0 and
max|E - E_profile| = 0.  Run this file directly for the self-check.
"""
import json
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / '02_人体骨骼晶格肌肉方案' / '输出' / 'anatomical_concept'
DEFINITION_PATH = OUTPUT_DIR / 'joint_model' / 'definition.json'
CANDIDATES_PATH = OUTPUT_DIR / 'joint_candidates.json'
MUSCLE_DIR = OUTPUT_DIR / 'elbow_muscles'

#: source -> display translation (display = source - SHIFT)
SHIFT = np.array([-140.0, -70.0, 1100.0])


def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError('zero-length vector')
    return v / n


def _load(path):
    with open(path, 'r', encoding='utf-8') as fh:
        return json.load(fh)


# ------------------------------------------------------------------ coordinates
def source_to_display(v):
    return np.asarray(v, dtype=float) - SHIFT


def display_to_source(v):
    return np.asarray(v, dtype=float) + SHIFT


def joint_centres(frame='source'):
    """Centres of SC / AC / GH / elbow.  frame: 'source' (default) or 'display'."""
    c = _load(DEFINITION_PATH)['centres_source_mm']
    out = {k: np.asarray(v, dtype=float) for k, v in c.items()}
    if frame == 'display':
        return {k: source_to_display(v) for k, v in out.items()}
    if frame != 'source':
        raise ValueError("frame must be 'source' or 'display'")
    return out


# ---------------------------------------------------------------------- frames
def muscle_frame():
    """Return (origin_display, basis_columns) of the muscle-local frame.

    `local = (point_display - origin) @ basis`.
    """
    c = joint_centres('display')
    E, G = c['elbow'], c['GH']
    z = unit(E - G)
    x = unit(np.cross(z, np.array([0.0, 1.0, 0.0])))
    y = np.cross(z, x)
    return E, np.stack([x, y, z], axis=1)


def gh_probe_axes():
    """The three GH axes used by audit_joint_model.py (rotation vectors, not Euler)."""
    c = joint_centres('source')
    long_axis = unit(c['elbow'] - c['GH'])
    return {
        'longitudinal': long_axis,                                   # axial rotation
        'flexion': unit(np.cross(long_axis, np.array([0.0, -1.0, 0.0]))),
        'abduction': unit(np.cross(long_axis, np.array([-1.0, 0.0, 0.0]))),
    }


def elbow_axis():
    """Return (centre_source, axis_source).  Positive angle = flexion, as probed in the project."""
    jc = _load(CANDIDATES_PATH)
    return (np.asarray(jc['elbow']['centre_source_mm'], dtype=float),
            unit(jc['elbow']['axis_source']))


# ------------------------------------------------------------------ self-check
def self_check(muscles=('FJ1486', 'FJ1477'), tol=1e-6):
    """Cross-check the derived frames against the stored muscle profiles."""
    report = {'pass': True, 'items': []}

    E, F = muscle_frame()
    orth = float(np.abs(F.T @ F - np.eye(3)).max())
    report['items'].append({'check': 'basis orthonormal', 'max|FtF-I|': orth,
                            'pass': orth < tol})
    report['pass'] &= orth < tol

    rt = float(np.abs(source_to_display(display_to_source(np.zeros(3)))).max())
    report['items'].append({'check': 'shift round-trip', 'error': rt, 'pass': rt < tol})
    report['pass'] &= rt < tol

    c = joint_centres('source')
    humerus = float(np.linalg.norm(c['elbow'] - c['GH']))
    report['items'].append({'check': 'GH-elbow humerus length mm', 'value': humerus})

    for mid in muscles:
        path = MUSCLE_DIR / (mid + '_profile.json')
        if not path.is_file():
            report['items'].append({'check': mid + ' profile', 'note': 'missing, skipped'})
            continue
        prof = _load(path)
        Ep = np.asarray(prof['coordinate_origin_display_mm'], dtype=float)
        Fp = np.asarray(prof['basis_columns_display'], dtype=float)
        de = float(np.abs(E - Ep).max())
        df = float(np.abs(F - Fp).max())
        ok = de < 1e-5 and df < 1e-5
        report['items'].append({'check': mid + ' profile frame', 'dE': de, 'dF': df, 'pass': ok})
        report['pass'] &= ok

    return report


def main():
    import sys
    rep = self_check()
    for item in rep['items']:
        flag = 'PASS' if item.get('pass', True) else 'FAIL'
        print('[%s] %s' % (flag, item))
    print('\nRESULT:', 'PASS' if rep['pass'] else 'FAIL')
    return 0 if rep['pass'] else 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
