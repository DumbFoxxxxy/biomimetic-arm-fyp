"""Turn 5.4.35's span changes into a requirement on the TPU: how much contraction is needed.

5.4.35 showed every axis has an agonist that shortens and an antagonist that lengthens, but it did
not ask whether the actuator CAN shorten that much.  The document records no achievable
contraction - that is exactly what S1-P0 and S1-G3 have to measure - so the useful thing to do is
invert it: state the requirement, and show how many muscles are needed if the bench measures less.

The required contraction of a line IS the relative span change of that line over the range (5.4.35),
so this reads that file and groups by axis.
"""
import json
import sys

import numpy as np

sys.path.insert(0, r'D:\XJTLU\PSP\mujoco_arm\tools')

import frames
import interference as I

CABLE = frames.OUTPUT_DIR.parent / 'shoulder_design' / 's1_cable_length.json'


def main():
    d = json.loads(CABLE.read_text(encoding='utf-8'))
    axes = ('flexion', 'abduction', 'longitudinal')
    summary = {}
    print('required contraction of the AGONISTS (a positive number), per axis:\n')
    for axis in axes:
        ag = sorted([r for r in d['rows'] if r['axis'] == axis and r['change_pct'] < 0],
                    key=lambda r: r['change_pct'])
        need = []
        for r in ag:
            need.append({'module': r['module'], 'head': r['head'],
                         'required_contraction_pct': round(-r['change_pct'], 2)})
        summary[axis] = need
        print('%-13s' % axis)
        for n in need:
            print('   %-14s %-8s needs %5.2f %%' % (n['module'], n['head'],
                                                    n['required_contraction_pct']))
        print()

    print('--- agonist trade-off: contraction vs moment arm ---')
    # A muscle with a tiny span change is cheap to contract but also produces almost no torque,
    # because d(span)/d(theta) IS the moment arm.  Both numbers have to be read together: dividing
    # the demand by the number of muscles (the first version of this script) says nothing, because
    # two muscles pulling the same DOF do not halve each other's stroke requirement.
    print('%-13s %-14s %-8s %12s %12s' %
          ('axis', 'module', 'head', 'contraction%', 'moment arm mm'))
    plan = {}
    for axis in axes:
        best = None
        for r in d['rows']:
            if r['axis'] != axis or r['change_pct'] >= 0:
                continue
            arm = abs(r['slope_mm_per_deg']) * 57.2958
            row = {'module': r['module'], 'head': r['head'],
                   'required_contraction_pct': round(-r['change_pct'], 2),
                   'moment_arm_mm': round(arm, 2)}
            plan.setdefault(axis, []).append(row)
            print('%-13s %-14s %-8s %12.2f %12.2f'
                  % (axis, r['module'], r['head'], -r['change_pct'], arm))
        plan[axis].sort(key=lambda x: x['required_contraction_pct'])

    # The axis requirement is a pair: a usable moment arm AND a feasible contraction.
    print('\n--- per-axis requirement (pick the agonist with arm >= 10 mm) ---')
    for axis in axes:
        usable = [r for r in plan[axis] if r['moment_arm_mm'] >= 10.0]
        if not usable:
            print('%-13s no agonist reaches a 10 mm moment arm' % axis)
            continue
        cheapest = min(usable, key=lambda r: r['required_contraction_pct'])
        print('%-13s %-14s %-8s contraction %.2f %%, arm %.1f mm'
              % (axis, cheapest['module'], cheapest['head'],
                 cheapest['required_contraction_pct'], cheapest['moment_arm_mm']))

    req = [{'axis': a, 'agonists': plan[a],
            'cheapest_usable_pct': (min([r['required_contraction_pct'] for r in plan[a]
                                         if r['moment_arm_mm'] >= 10.0])
                                    if any(r['moment_arm_mm'] >= 10.0 for r in plan[a]) else None)}
           for a in axes]
    print('\nreport scope: the lattice must achieve at least the cheapest usable agonist above;')
    print('the actual achievable contraction is what S1-P0 and S1-G3 must measure.')
    I.write_report(str(frames.OUTPUT_DIR.parent / 'shoulder_design' /
                       's1_contraction_requirement.json'), {
        'scope': 'per-agonist required contraction for the full 0-120 deg plus its moment arm '
                 '(= d(span)/d(theta)); the achievable contraction is what S1-P0/S1-G3 measure',
        'axes': req,
        'note': 'contraction and moment arm trade off: a muscle whose span barely changes is cheap '
                'to contract and useless as a torque source. Do NOT divide the demand by the '
                'number of muscles - two muscles on one DOF do not halve each stroke.',
        'not_verified': ['no achievable contraction is known yet - these are REQUIREMENTS',
                         'the share-by-k idea was wrong and is recorded here as retracted',
                         'chord span, not corridor length',
                         'no series compliance, no tendon stiffness'],
    })
    print('\nreport written: s1_contraction_requirement.json')


if __name__ == '__main__':
    main()
